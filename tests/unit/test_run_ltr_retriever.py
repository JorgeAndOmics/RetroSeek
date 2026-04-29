"""Unit tests for ``workflow/scripts/run_ltr_retriever.py``.

The runner wraps LTR_retriever's Perl pipeline so the Snakemake rule
becomes a one-line `python` call. Tests cover the four building-block
helpers (resolve_source_scn, stage_workdir, run_binary, finalise_outputs)
plus an end-to-end main() flow with a bash shim binary.
"""

from __future__ import annotations

import stat
from pathlib import Path

import pytest

from run_ltr_retriever import (
    EXPECTED_OUTPUT_EXTS,
    LTRRetrieverParams,
    finalise_outputs,
    main,
    resolve_source_scn,
    run_binary,
    stage_workdir,
)


# ---------------------------------------------------------------------
# resolve_source_scn
# ---------------------------------------------------------------------
def test_resolve_source_scn_picks_retroviral_when_configured(tmp_path: Path) -> None:
    retroviral = tmp_path / "r.scn"
    full = tmp_path / "f.scn"
    retroviral.write_text("retroviral\n")
    full.write_text("full\n")
    assert resolve_source_scn("retroviral", retroviral, full) == retroviral


def test_resolve_source_scn_picks_full_when_configured(tmp_path: Path) -> None:
    retroviral = tmp_path / "r.scn"
    full = tmp_path / "f.scn"
    retroviral.write_text("retroviral\n")
    full.write_text("full\n")
    assert resolve_source_scn("full", retroviral, full) == full


def test_resolve_source_scn_raises_on_unknown_mode(tmp_path: Path) -> None:
    """Bad config should fail loudly, not default silently."""
    with pytest.raises(ValueError, match="source_scn"):
        resolve_source_scn("banana", tmp_path / "r.scn", tmp_path / "f.scn")


# ---------------------------------------------------------------------
# stage_workdir
# ---------------------------------------------------------------------
def _make_genome_and_scn(tmp_path: Path) -> tuple[Path, Path]:
    genome = tmp_path / "src" / "Toyus.fa"
    genome.parent.mkdir(parents=True)
    genome.write_text(">chr1\nACGT\n")
    scn = tmp_path / "src" / "Toyus_retroviral.scn"
    scn.write_text("# header\n100 200 ...\n")
    return genome, scn


def test_stage_workdir_creates_clean_dir_with_symlinks(tmp_path: Path) -> None:
    genome, scn = _make_genome_and_scn(tmp_path)
    workdir = tmp_path / "workdir"

    staged = stage_workdir(workdir, genome, scn, genome_name="Toyus")

    assert workdir.is_dir()
    assert staged.workdir == workdir
    assert staged.genome_name == "Toyus"
    fa_link = workdir / "Toyus.fa"
    scn_link = workdir / "Toyus.scn"
    assert fa_link.is_symlink()
    assert fa_link.resolve() == genome.resolve()
    assert scn_link.is_symlink()
    assert scn_link.resolve() == scn.resolve()


def test_stage_workdir_replaces_stale_symlink(tmp_path: Path) -> None:
    """A pre-existing symlink to the wrong target must be replaced cleanly."""
    genome, scn = _make_genome_and_scn(tmp_path)
    workdir = tmp_path / "workdir"
    workdir.mkdir()
    # Plant a stale symlink pointing to a nonexistent path.
    stale = workdir / "Toyus.fa"
    stale.symlink_to(tmp_path / "does_not_exist")
    assert stale.is_symlink()

    stage_workdir(workdir, genome, scn, genome_name="Toyus")

    assert (workdir / "Toyus.fa").resolve() == genome.resolve()


# ---------------------------------------------------------------------
# run_binary
# ---------------------------------------------------------------------
def _basic_params(binary: Path = Path("/bin/echo")) -> LTRRetrieverParams:
    return LTRRetrieverParams(
        substitution_rate=1.3e-8,
        min_similarity=91,
        threads=1,
        noanno=True,
        binary=binary,
    )


def test_run_binary_captures_stdout_and_stderr_to_log(tmp_path: Path) -> None:
    """stdout + stderr are tee'd to the log file with prefixes."""
    genome, scn = _make_genome_and_scn(tmp_path)
    workdir = tmp_path / "wd"
    staged = stage_workdir(workdir, genome, scn, genome_name="Toyus")

    log = tmp_path / "run.log"

    # Use a tiny shell shim that writes one line each to stdout/stderr.
    shim = tmp_path / "shim.sh"
    shim.write_text(
        "#!/bin/bash\necho 'hello on stdout'\necho 'hello on stderr' >&2\nexit 0\n"
    )
    shim.chmod(shim.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    rc = run_binary(staged, _basic_params(binary=shim), log)
    assert rc == 0
    text = log.read_text()
    assert "hello on stdout" in text
    assert "hello on stderr" in text


def test_run_binary_propagates_nonzero_exit_code(tmp_path: Path) -> None:
    genome, scn = _make_genome_and_scn(tmp_path)
    workdir = tmp_path / "wd"
    staged = stage_workdir(workdir, genome, scn, genome_name="Toyus")

    shim = tmp_path / "fail.sh"
    shim.write_text("#!/bin/bash\necho 'about to fail' >&2\nexit 7\n")
    shim.chmod(shim.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    log = tmp_path / "fail.log"
    rc = run_binary(staged, _basic_params(binary=shim), log)
    assert rc == 7
    # Log must still be written even on failure.
    assert "about to fail" in log.read_text()


# ---------------------------------------------------------------------
# finalise_outputs
# ---------------------------------------------------------------------
def test_finalise_outputs_renames_fa_mod_prefix(tmp_path: Path) -> None:
    genome_name = "Toyus"
    for ext in EXPECTED_OUTPUT_EXTS:
        (tmp_path / f"{genome_name}.fa.mod.{ext}").write_text(f"data {ext}\n")

    canonical = finalise_outputs(tmp_path, genome_name)

    assert sorted(p.name for p in canonical) == sorted(
        f"{genome_name}.{ext}" for ext in EXPECTED_OUTPUT_EXTS
    )
    for ext in EXPECTED_OUTPUT_EXTS:
        assert (tmp_path / f"{genome_name}.{ext}").is_file()
        assert not (tmp_path / f"{genome_name}.fa.mod.{ext}").exists()


def test_finalise_outputs_raises_on_missing_expected_file(tmp_path: Path) -> None:
    """Only 2 of 3 expected outputs present → loud failure with a useful message."""
    genome_name = "Toyus"
    # Materialise only two of the three.
    (tmp_path / f"{genome_name}.fa.mod.pass.list.gff3").write_text("x")
    (tmp_path / f"{genome_name}.fa.mod.LTRlib.fa").write_text("x")

    with pytest.raises(RuntimeError, match="nmtf.pass.list"):
        finalise_outputs(tmp_path, genome_name)


def test_finalise_outputs_handles_already_canonical_files(tmp_path: Path) -> None:
    """If the canonical names already exist (no .fa.mod prefix), idempotent."""
    genome_name = "Toyus"
    for ext in EXPECTED_OUTPUT_EXTS:
        (tmp_path / f"{genome_name}.{ext}").write_text(f"already {ext}\n")

    canonical = finalise_outputs(tmp_path, genome_name)
    assert len(canonical) == len(EXPECTED_OUTPUT_EXTS)


# ---------------------------------------------------------------------
# main — end-to-end with a fake LTR_retriever binary
# ---------------------------------------------------------------------
@pytest.mark.integration
def test_main_end_to_end_with_fake_binary(tmp_path: Path) -> None:
    """Full main() flow against a bash shim that emits .fa.mod.* outputs."""
    genome, scn_full = _make_genome_and_scn(tmp_path)
    scn_retroviral = tmp_path / "src" / "Toyus_retroviral.scn"
    scn_retroviral.write_text("# header\n100 200 ...\n")

    workdir = tmp_path / "wd"
    log = tmp_path / "main.log"

    shim = tmp_path / "fake_LTR_retriever.sh"
    shim.write_text(
        "#!/bin/bash\n"
        "echo 'fake LTR_retriever running'\n"
        "echo 'genome flag was:' $2\n"
        # LTR_retriever's working dir is the cwd; produce the .fa.mod.* trio.
        "for ext in pass.list.gff3 nmtf.pass.list LTRlib.fa; do\n"
        "  echo 'fake content' > Toyus.fa.mod.$ext\n"
        "done\n"
        "exit 0\n"
    )
    shim.chmod(shim.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    rc = main(
        [
            "--genome-fa", str(genome),
            "--retroviral-scn", str(scn_retroviral),
            "--full-scn", str(scn_full),
            "--source-scn-mode", "retroviral",
            "--workdir", str(workdir),
            "--genome-name", "Toyus",
            "--substitution-rate", "1.3e-8",
            "--min-similarity", "91",
            "--threads", "1",
            "--noanno",
            "--log-file", str(log),
            "--ltr-retriever-binary", str(shim),
        ]
    )

    assert rc == 0
    for ext in EXPECTED_OUTPUT_EXTS:
        assert (workdir / f"Toyus.{ext}").is_file()
    assert "fake LTR_retriever running" in log.read_text()


def test_main_returns_nonzero_when_binary_fails(tmp_path: Path) -> None:
    """Failing LTR_retriever exits non-zero from main; no output renaming attempted."""
    genome, scn_full = _make_genome_and_scn(tmp_path)
    workdir = tmp_path / "wd"
    log = tmp_path / "main.log"

    shim = tmp_path / "fail.sh"
    shim.write_text("#!/bin/bash\necho 'boom' >&2\nexit 3\n")
    shim.chmod(shim.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    rc = main(
        [
            "--genome-fa", str(genome),
            "--retroviral-scn", str(scn_full),
            "--full-scn", str(scn_full),
            "--source-scn-mode", "retroviral",
            "--workdir", str(workdir),
            "--genome-name", "Toyus",
            "--substitution-rate", "1.3e-8",
            "--min-similarity", "91",
            "--threads", "1",
            "--log-file", str(log),
            "--ltr-retriever-binary", str(shim),
        ]
    )
    assert rc == 3
