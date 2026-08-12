"""Run LTR_retriever in a staged workdir, capture logs, normalise outputs.

Why this script exists
----------------------
LTR_retriever is a Perl pipeline that:

1. expects its working directory to be the directory containing the
   genome FASTA - it writes outputs alongside the input;
2. produces files prefixed with ``<basename>.fa.mod.`` because it
   internally generates a sanitised copy of the genome before running;
3. emits no machine-readable success signal - the Snakemake rule has
   to verify the three expected output files materialised.

Embedding all that in a Snakemake ``shell:`` block produced fragile
multi-line bash that silently no-op'd on missing files. This wrapper
isolates each concern into a unit-testable Python helper.

Where solo LTRs actually come from
----------------------------------
LTR_retriever does **not** report solo LTRs in ``nmtf.pass.list``. That file
holds *intact* LTR-RTs whose termini lack the canonical TGCA motif - the
tool's own result banner calls it ``(Non-TGCA LTR-RTs)`` and its summary line
reads "Total intact non-TGCA LTR-RTs found". An earlier revision of this
workstream read "nmtf" as "non-matching-full" and wired the integrator to it,
which would have reported intact elements as solo LTRs.

The real path runs off the whole-genome RepeatMasker annotation:

1. LTR_retriever annotates the genome with its own LTR library, writing
   ``{genome}.out`` (RepeatMasker table). This only happens when annotation is
   enabled, so this runner never passes ``-noanno``.
2. ``bin/find_LTR.pl -lib {genome}.LTRlib.fa`` maps the LTR regions inside each
   library sequence.
3. ``bin/solo_finder.pl -i {genome}.out -info {genome}.LTR.info`` emits the solo
   list: ``chrom, start, end, locus, library_id, coverage``. A hit counts as
   solo when it covers 0.8-1.2 of the library LTR, is at least 80 bp, and sits
   at least 300 bp clear of any internal-region annotation.

CLI
---
::

    python run_ltr_retriever.py \\
        --genome-fa <path> \\
        --retroviral-scn <path> --full-scn <path> \\
        --source-scn-mode <retroviral|full> \\
        --workdir <path> \\
        --genome-name <str> \\
        --substitution-rate <float> --min-similarity <int> --threads <int> \\
        --log-file <path> [--ltr-retriever-binary <path>]
"""

from __future__ import annotations

import argparse
import logging
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

# `{genome}.out` is RepeatMasker's whole-genome table and the sole input to
# solo_finder.pl; `pass.list` is the intact set that forms the solo/intact
# ratio's denominator. Both are absent unless annotation runs.
EXPECTED_OUTPUT_EXTS: tuple[str, ...] = (
    "pass.list",
    "pass.list.gff3",
    "LTRlib.fa",
    "out",
)
VALID_SOURCE_SCN_MODES = ("retroviral", "full")


@dataclass
class StagedPaths:
    """Resolved paths inside the per-genome workdir."""

    workdir: Path
    genome_fa_link: Path
    scn_link: Path
    genome_name: str


@dataclass
class LTRRetrieverParams:
    """Tunable parameters forwarded to the LTR_retriever binary.

    ``-noanno`` is deliberately absent: it suppresses the whole-genome
    RepeatMasker pass, and that pass produces the only file solo_finder.pl can
    read. Solo-LTR detection is the entire purpose of this stage, so annotation
    is not optional here.
    """

    substitution_rate: float
    min_similarity: int
    threads: int
    binary: Path


# ---------------------------------------------------------------------
# resolve_source_scn
# ---------------------------------------------------------------------
def resolve_source_scn(mode: str, retroviral: Path, full: Path) -> Path:
    """Pick the SCN that should feed LTR_retriever.

    ``mode`` mirrors ``config.ltr_retriever.source_scn``:

    - ``retroviral`` - the prefilter-restricted SCN (Coupling A).
    - ``full`` - the unfiltered passthrough SCN.

    Any other value raises ``ValueError`` rather than silently
    defaulting; the validator should have caught it but we double-
    check at the runner boundary.
    """
    if mode == "retroviral":
        return retroviral
    if mode == "full":
        return full
    raise ValueError(
        f"Unknown source_scn mode {mode!r}; expected one of {VALID_SOURCE_SCN_MODES}"
    )


# ---------------------------------------------------------------------
# stage_workdir
# ---------------------------------------------------------------------
def _ensure_symlink(link: Path, target: Path) -> None:
    """Make ``link`` point at ``target``, replacing any stale link in place."""
    if link.is_symlink() or link.exists():
        link.unlink()
    link.symlink_to(target.resolve())


def stage_workdir(
    workdir: Path,
    genome_fa: Path,
    scn: Path,
    genome_name: str,
) -> StagedPaths:
    """Materialise a clean per-genome workdir with symlinks to inputs.

    LTR_retriever runs in this directory and writes its outputs here.
    Symlinks (not copies) keep gigabyte-scale genomes off the
    intermediate filesystem; LTR_retriever doesn't modify the input.
    Stale symlinks from prior runs are replaced rather than reused.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    fa_link = workdir / f"{genome_name}.fa"
    scn_link = workdir / f"{genome_name}.scn"
    _ensure_symlink(fa_link, genome_fa)
    _ensure_symlink(scn_link, scn)
    return StagedPaths(
        workdir=workdir,
        genome_fa_link=fa_link,
        scn_link=scn_link,
        genome_name=genome_name,
    )


# ---------------------------------------------------------------------
# run_binary
# ---------------------------------------------------------------------
def _build_command(
    staged: StagedPaths,
    params: LTRRetrieverParams,
) -> list[str]:
    """Construct the LTR_retriever argv list."""
    return [
        str(params.binary),
        "-genome",
        staged.genome_fa_link.name,
        "-inharvest",
        staged.scn_link.name,
        "-u",
        str(params.substitution_rate),
        "-miniden",
        str(params.min_similarity),
        "-threads",
        str(params.threads),
    ]


def run_binary(
    staged: StagedPaths,
    params: LTRRetrieverParams,
    log_file: Path,
) -> int:
    """Invoke LTR_retriever, tee stdout+stderr to ``log_file``.

    Returns the binary's exit code. On nonzero exit, the log is left
    intact so the caller (or the user) can inspect it.
    """
    log_file.parent.mkdir(parents=True, exist_ok=True)
    cmd = _build_command(staged, params)
    with log_file.open("w") as log:
        log.write(f"# LTR_retriever runner\n# cmd: {' '.join(cmd)}\n")
        log.write(f"# cwd: {staged.workdir}\n\n")
        log.flush()
        proc = subprocess.run(
            cmd,
            cwd=staged.workdir,
            capture_output=True,
            check=False,
            text=True,
        )
        for line in proc.stdout.splitlines():
            log.write(f"[stdout] {line}\n")
        for line in proc.stderr.splitlines():
            log.write(f"[stderr] {line}\n")
        log.write(f"\n# exit code: {proc.returncode}\n")
    return proc.returncode


# ---------------------------------------------------------------------
# finalise_outputs
# ---------------------------------------------------------------------
def finalise_outputs(workdir: Path, genome_name: str) -> list[Path]:
    """Rename ``{genome}.fa.mod.<ext>`` -> ``{genome}.<ext>`` for each expected ext.

    If the canonical filename already exists (no ``.fa.mod.`` prefix),
    it is left in place. If neither the prefixed nor the canonical
    file exists for any expected extension, raise ``RuntimeError``
    with a list of missing files and a workdir directory listing -
    the strongest signal that LTR_retriever failed silently.
    """
    canonical: list[Path] = []
    missing: list[str] = []
    for ext in EXPECTED_OUTPUT_EXTS:
        canonical_path = workdir / f"{genome_name}.{ext}"
        prefixed_path = workdir / f"{genome_name}.fa.mod.{ext}"
        if canonical_path.exists():
            canonical.append(canonical_path)
            # If the prefixed sibling also exists, drop the redundant copy.
            if prefixed_path.exists():
                prefixed_path.unlink()
        elif prefixed_path.exists():
            prefixed_path.rename(canonical_path)
            canonical.append(canonical_path)
        else:
            missing.append(f"{genome_name}.{ext}")
    if missing:
        listing = sorted(p.name for p in workdir.iterdir())
        raise RuntimeError(
            f"LTR_retriever finished but expected outputs missing: {missing}. "
            f"Workdir contents: {listing}"
        )
    return canonical


# ---------------------------------------------------------------------
# solo finding
# ---------------------------------------------------------------------
def resolve_helper_dir(binary: Path) -> Path:
    """Locate LTR_retriever's ``bin/`` directory of Perl helper scripts.

    ``find_LTR.pl`` and ``solo_finder.pl`` are not installed on PATH; they live
    beside the main Perl program under ``share/LTR_retriever/bin``. Conda ships
    ``bin/LTR_retriever`` as a two-line bash shim that execs the real script, so
    the shim's own directory has no ``bin/`` subdirectory. Read the shim to find
    the interpreter target when it is not the Perl program itself.

    Raises
    ------
    FileNotFoundError
        If no directory containing the helper scripts can be found.
    """
    candidates = [binary.parent / "bin"]
    try:
        text = binary.read_text(errors="replace")
    except OSError:
        text = ""
    candidates.extend(
        Path(match.rstrip("$@ ")).parent / "bin"
        for match in re.findall(r"(\S*share/LTR_retriever\S*)", text)
    )
    candidates.append(binary.parent.parent / "share" / "LTR_retriever" / "bin")
    for candidate in candidates:
        if (candidate / "solo_finder.pl").is_file():
            return candidate
    raise FileNotFoundError(
        "Could not locate LTR_retriever's bin/ helper scripts (solo_finder.pl) "
        f"from binary {binary}. Looked in: {[str(c) for c in candidates]}"
    )


def run_solo_finder(
    staged: StagedPaths,
    helper_dir: Path,
    log_file: Path,
) -> Path:
    """Derive the solo-LTR list from LTR_retriever's whole-genome annotation.

    Chains the tool's own two helpers rather than reimplementing their criteria:
    ``find_LTR.pl`` reports where the LTR regions sit inside each library
    sequence, and ``solo_finder.pl`` walks the RepeatMasker table keeping hits
    that look like a lone LTR rather than one flank of an intact element.

    Returns the path to ``{genome}.solo_list``. An empty file is a legitimate
    result for a genome with no solos, so it is not treated as failure - but a
    missing ``{genome}.out`` is, since that means annotation never ran.
    """
    genome = staged.genome_name
    ltr_lib = staged.workdir / f"{genome}.LTRlib.fa"
    rm_out = staged.workdir / f"{genome}.out"
    ltr_info = staged.workdir / f"{genome}.LTR.info"
    solo_list = staged.workdir / f"{genome}.solo_list"

    if not rm_out.is_file():
        raise RuntimeError(
            f"RepeatMasker output {rm_out} is missing, so solo LTRs cannot be "
            "called. LTR_retriever must run with whole-genome annotation "
            "enabled (this runner never passes -noanno)."
        )

    steps = (
        (["perl", str(helper_dir / "find_LTR.pl"), "-lib", str(ltr_lib)], ltr_info),
        (
            [
                "perl",
                str(helper_dir / "solo_finder.pl"),
                "-i",
                str(rm_out),
                "-info",
                str(ltr_info),
            ],
            solo_list,
        ),
    )
    with log_file.open("a") as log:
        for cmd, destination in steps:
            log.write(f"\n# solo step: {' '.join(cmd)} > {destination.name}\n")
            log.flush()
            proc = subprocess.run(
                cmd, cwd=staged.workdir, capture_output=True, check=False, text=True
            )
            for line in proc.stderr.splitlines():
                log.write(f"[stderr] {line}\n")
            if proc.returncode != 0:
                log.write(f"# exit code: {proc.returncode}\n")
                raise RuntimeError(f"{cmd[1]} exited {proc.returncode}; see {log_file}")
            destination.write_text(proc.stdout)
    return solo_list


# ---------------------------------------------------------------------
# main
# ---------------------------------------------------------------------
def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--genome-fa", type=Path, required=True)
    parser.add_argument("--retroviral-scn", type=Path, required=True)
    parser.add_argument("--full-scn", type=Path, required=True)
    parser.add_argument(
        "--source-scn-mode",
        choices=VALID_SOURCE_SCN_MODES,
        required=True,
    )
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--genome-name", required=True)
    parser.add_argument("--substitution-rate", type=float, required=True)
    parser.add_argument("--min-similarity", type=int, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--log-file", type=Path, required=True)
    parser.add_argument("--ltr-retriever-binary", type=Path, default=None)
    return parser.parse_args(argv)


def _resolve_binary(explicit: Path | None) -> Path:
    if explicit is not None:
        return explicit
    found = shutil.which("LTR_retriever")
    if found is None:
        raise FileNotFoundError(
            "LTR_retriever not found on PATH and --ltr-retriever-binary not given"
        )
    return Path(found)


def main(argv: list[str] | None = None) -> int:
    """Entry point - orchestrates the four building blocks."""
    args = _parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    scn = resolve_source_scn(args.source_scn_mode, args.retroviral_scn, args.full_scn)
    if not args.genome_fa.exists():
        raise FileNotFoundError(f"genome FASTA not found: {args.genome_fa}")
    if not scn.exists():
        raise FileNotFoundError(f"selected SCN not found: {scn}")

    binary = _resolve_binary(args.ltr_retriever_binary)
    params = LTRRetrieverParams(
        substitution_rate=args.substitution_rate,
        min_similarity=args.min_similarity,
        threads=args.threads,
        binary=binary,
    )

    staged = stage_workdir(args.workdir, args.genome_fa, scn, args.genome_name)
    logger.info(
        "LTR_retriever: genome=%s scn=%s mode=%s workdir=%s",
        args.genome_fa.name,
        scn.name,
        args.source_scn_mode,
        args.workdir,
    )
    rc = run_binary(staged, params, args.log_file)
    if rc != 0:
        logger.error("LTR_retriever exited %d; see %s", rc, args.log_file)
        return rc

    finalise_outputs(staged.workdir, args.genome_name)
    solo_list = run_solo_finder(staged, resolve_helper_dir(binary), args.log_file)
    n_solo = sum(1 for line in solo_list.read_text().splitlines() if line.strip())
    logger.info("solo_finder reported %d solo LTR(s) -> %s", n_solo, solo_list)
    return 0


if __name__ == "__main__":
    sys.exit(main())
