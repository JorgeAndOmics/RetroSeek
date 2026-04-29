"""Run LTR_retriever in a staged workdir, capture logs, normalise outputs.

Why this script exists
----------------------
LTR_retriever is a Perl pipeline that:

1. expects its working directory to be the directory containing the
   genome FASTA — it writes outputs alongside the input;
2. produces files prefixed with ``<basename>.fa.mod.`` because it
   internally generates a sanitised copy of the genome before running;
3. emits no machine-readable success signal — the Snakemake rule has
   to verify the three expected output files materialised.

Embedding all that in a Snakemake ``shell:`` block produced fragile
multi-line bash that silently no-op'd on missing files. This wrapper
isolates each concern into a unit-testable Python helper.

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
        [--noanno] \\
        --log-file <path> [--ltr-retriever-binary <path>]
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

EXPECTED_OUTPUT_EXTS: tuple[str, ...] = (
    "pass.list.gff3",
    "nmtf.pass.list",
    "LTRlib.fa",
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
    """Tunable parameters forwarded to the LTR_retriever binary."""

    substitution_rate: float
    min_similarity: int
    threads: int
    noanno: bool
    binary: Path


# ---------------------------------------------------------------------
# resolve_source_scn
# ---------------------------------------------------------------------
def resolve_source_scn(mode: str, retroviral: Path, full: Path) -> Path:
    """Pick the SCN that should feed LTR_retriever.

    ``mode`` mirrors ``config.ltr_retriever.source_scn``:

    - ``retroviral`` — the prefilter-restricted SCN (Coupling A).
    - ``full`` — the unfiltered passthrough SCN.

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
    cmd = [
        str(params.binary),
        "-genome", staged.genome_fa_link.name,
        "-inharvest", staged.scn_link.name,
        "-u", str(params.substitution_rate),
        "-miniden", str(params.min_similarity),
        "-threads", str(params.threads),
    ]
    if params.noanno:
        cmd.append("-noanno")
    return cmd


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
        proc = subprocess.run(  # noqa: S603 — args are constructed from validated paths
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
    """Rename ``{genome}.fa.mod.<ext>`` → ``{genome}.<ext>`` for each expected ext.

    If the canonical filename already exists (no ``.fa.mod.`` prefix),
    it is left in place. If neither the prefixed nor the canonical
    file exists for any expected extension, raise ``RuntimeError``
    with a list of missing files and a workdir directory listing —
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
    parser.add_argument("--noanno", action="store_true")
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
    """Entry point — orchestrates the four building blocks."""
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
        noanno=args.noanno,
        binary=binary,
    )

    staged = stage_workdir(args.workdir, args.genome_fa, scn, args.genome_name)
    logging.info(
        "LTR_retriever: genome=%s scn=%s mode=%s workdir=%s",
        args.genome_fa.name,
        scn.name,
        args.source_scn_mode,
        args.workdir,
    )
    rc = run_binary(staged, params, args.log_file)
    if rc != 0:
        logging.error("LTR_retriever exited %d; see %s", rc, args.log_file)
        return rc

    finalise_outputs(staged.workdir, args.genome_name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
