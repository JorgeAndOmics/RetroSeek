# =============================================================================
# subset_pfam.py
# =============================================================================
# Derive a small, curated HMM library from the full Pfam-A flat file.
#
# Why this exists: the domain scan needs one library applied identically to
# LTR-flanked and orphan loci. Searching all 30,134 Pfam models over every locus
# costs hours; searching the ~150 curated families costs minutes. At Pfam's
# gathering thresholds (`hmmsearch --cut_ga`) only 39 distinct families fire
# across the orphan marker sets, so the subset is very nearly lossless - the
# 6,175-family "tail" LTRdigest reports is an artifact of running with no
# per-family cutoff.
#
# Two invariants, both load-bearing:
#   1. Selection is by ACCESSION, never by name. Pfam guarantees accession
#      stability (PF00665 is integrase forever) but reserves the right to rename
#      families. The file stores versioned accessions (`ACC PF00665.33`); the
#      curated table stores the unversioned key, so the version is stripped here.
#   2. An accession the table asks for and the library does not contain is a
#      hard error naming the offenders, the likely cause and the fix. Silently
#      emitting a short library would make every downstream count wrong with no
#      signal. The usual cause is a Pfam-A.hmm older than the pinned release
#      (ADR-019): new families such as PF29843 do not exist in it yet.
#
# The file is read as a stream (one record at a time), not slurped: Pfam-A.hmm is
# ~2.2 GB. Measured on the full file: 5.3 s, 116 MB peak RSS.
#
# `write_name_map` is separate and reads the FULL library, because it exists to
# resolve the domain NAMEs that `gt ltrdigest` writes (it never emits accessions)
# back to accessions, including families the curated table does not list.
# =============================================================================

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from log import OK, PipelineError, job_logging, run_main

logger = logging.getLogger(__name__)


def wanted_accessions(classes_tsv: Path) -> set[str]:
    """Unversioned Pfam accessions named by the curated class table.

    The table is `pfam_acc <TAB> pfam_name <TAB> class` with a header row.
    """
    wanted: set[str] = set()
    with classes_tsv.open(encoding="utf-8") as fh:
        next(fh, None)  # header
        for raw in fh:
            row = raw.strip()
            if not row or row.startswith("#"):
                continue
            wanted.add(row.split("\t")[0].split(".")[0])
    return wanted


def library_accessions(hmm_path: Path) -> set[str]:
    """Unversioned accessions of every model in `hmm_path`, read as a stream."""
    found: set[str] = set()
    with hmm_path.open(encoding="utf-8") as source:
        for line in source:
            if line.startswith("ACC "):
                found.add(line.split()[1].split(".")[0])
    return found


def missing_accessions(hmm_path: Path, wanted: set[str]) -> list[str]:
    """The accessions in `wanted` that `hmm_path` does not hold, sorted.

    Used before a run starts (the launcher's preflight), so a stale library is
    reported in seconds instead of after the upstream stages have run.
    """
    return sorted(wanted - library_accessions(hmm_path))


MISSING_HINT = (
    "./RetroSeek --download-hmm --forcerun pfam_hmm_downloader, then backdate the "
    "new files so LTRdigest does not rerun (docs/adr/ADR-019)"
)


def missing_message(missing: list[str], hmm_path: Path) -> str:
    """Explain a library that lacks curated accessions: what, why, and the fix."""
    return (
        f"{hmm_path} lacks {len(missing)} accession(s) listed in the curated "
        f"Pfam class table: {', '.join(missing)}. The file probably comes from an "
        "older Pfam release than the one the pipeline pins (input.pfam_release). "
        f"Fix: {MISSING_HINT}."
    )


def subset_pfam(hmm_path: Path, classes_tsv: Path, out_path: Path) -> int:
    """Copy the records named by `classes_tsv` from `hmm_path` into `out_path`.

    Returns the number of models written. Exits naming any accession that the
    table requested and the library does not hold.
    """
    wanted = wanted_accessions(classes_tsv)
    if not wanted:
        raise PipelineError(
            f"no accessions listed in {classes_tsv}",
            hint="point input.pfam_domain_classes at the curated class table",
        )

    found: set[str] = set()
    record: list[str] = []
    accession: str | None = None
    with (
        hmm_path.open(encoding="utf-8") as source,
        out_path.open("w", encoding="utf-8") as out,
    ):
        for line in source:
            record.append(line)
            if line.startswith("ACC "):
                accession = line.split()[1].split(".")[0]
            elif line.startswith("//"):
                if accession in wanted:
                    out.writelines(record)
                    found.add(accession)
                record, accession = [], None

    missing = sorted(wanted - found)
    if missing:
        raise PipelineError(
            f"{hmm_path} lacks {len(missing)} accession(s) listed in the curated "
            f"Pfam class table: {', '.join(missing)}. The file probably comes from "
            "an older Pfam release than the one the pipeline pins "
            "(input.pfam_release)",
            hint=MISSING_HINT,
        )
    return len(found)


def write_name_map(hmm_path: Path, out_path: Path) -> int:
    """Write a `pfam_name <TAB> pfam_acc` map covering every model in `hmm_path`.

    Needed because LTRdigest records only the model NAME. Names are unique within
    one Pfam release (38.2: 30,134 names, 30,134 accessions, no duplicates) but
    not across releases, so this map is only valid for the library it was built
    from - which is why it is generated beside the subset rather than committed.
    """
    written = 0
    name: str | None = None
    with (
        hmm_path.open(encoding="utf-8") as source,
        out_path.open("w", encoding="utf-8") as out,
    ):
        for line in source:
            if line.startswith("NAME "):
                name = line.split(maxsplit=1)[1].strip()
            elif line.startswith("ACC ") and name is not None:
                out.write(f"{name}\t{line.split()[1].split('.')[0]}\n")
                written += 1
                name = None
    return written


def main(args: argparse.Namespace) -> None:
    n_models = subset_pfam(args.hmm, args.classes, args.out)
    n_names = write_name_map(args.hmm, args.name_map)
    logger.info("wrote %d name to accession pairs to %s", n_names, args.name_map)
    logger.log(OK, "%s curated Pfam models in the subset", f"{n_models:,}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build the curated Pfam subset.")
    parser.add_argument("--hmm", type=Path, required=True, help="full Pfam-A.hmm")
    parser.add_argument("--classes", type=Path, required=True, help="curated class TSV")
    parser.add_argument("--out", type=Path, required=True, help="subset HMM to write")
    parser.add_argument(
        "--name-map", type=Path, required=True, help="name to accession map to write"
    )
    parser.add_argument("--log", type=Path, help="job log (the Snakemake log: path)")
    cli = parser.parse_args()
    job_logging(cli.log, "pfam_subset_builder")
    run_main(lambda: main(cli))
