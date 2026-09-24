"""Build the solo-LTR search bait: the LTR arms of ERV-bearing elements.

This is the step that makes the rest of the method classifier-free. An LTR arm
belonging to an element that hosts a catalogued ERV locus is retroviral *by
construction*: the probe search already found retroviral protein inside that
element, so its LTRs are the LTRs of a retrovirus. Using those arms as blastn bait
means every hit is a copy of a known-retroviral LTR, and nothing downstream has to
decide whether a sequence "looks retroviral".

Inputs:
    the flanking-LTR track (`tracks/flanking_ltr/{genome}.gff3`), which holds every
    LTRharvest arm as a `long_terminal_repeat` with `Parent=` and `arm=`;
    the classified loci table (`{genome}.loci.csv`), whose `parent` column names the
    element hosting each catalogued LTR-flanked ERV locus.

Output:
    a BED6 of the selected arms, ready for `taxonomy/extract_region_fasta.R`.

Arms below `--min-bait-length` are dropped rather than trusted. The coverage rule
that accepts a hit is a *fraction* of the bait, and 80% of a 102 bp arm is 82 bp,
which in a genome several percent LTR by mass is not evidence of anything. That one
omission is what produced the prototype's 492:1 solo/intact ratio (ADR-017).
"""

from __future__ import annotations

import argparse
import csv
import logging
from collections.abc import Iterator
from pathlib import Path
from typing import NamedTuple

from log import OK, PipelineError, job_logging, run_main

logger = logging.getLogger(__name__)

ARM_FEATURE = "long_terminal_repeat"


class Arm(NamedTuple):
    """One LTR arm, in GFF3 coordinates (1-based, both ends inclusive)."""

    seqname: str
    start: int
    end: int
    parent: str
    arm: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    @property
    def name(self) -> str:
        """The FASTA/BED name, and the only route back from a hit to its element."""
        return f"{self.seqname}|{self.parent}|{self.arm}"


def _attributes(column: str) -> dict[str, str]:
    """Parse a GFF3 attributes column into a dict, ignoring malformed entries."""
    out = {}
    for field in column.rstrip().split(";"):
        key, _, value = field.partition("=")
        if value:
            out[key.strip()] = value.strip()
    return out


def erv_bearing_parents(loci_csv: Path) -> set[str]:
    """Element IDs that host a catalogued ERV locus, from the loci table.

    Only LTR-flanked loci have a parent element, so no source filter is needed: the
    orphan tier lives in its own file. A missing `parent` column is an error rather
    than an empty result, because "no bait" and "no ERVs in this genome" would
    otherwise be indistinguishable.
    """
    with loci_csv.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "parent" not in reader.fieldnames:
            raise PipelineError(
                f"{loci_csv} has no 'parent' column, so no element can be identified "
                f"as ERV-bearing (columns present: {reader.fieldnames})",
                hint="rerun --classify; the loci table predates the parent column",
            )
        return {row["parent"] for row in reader if row.get("parent")}


def parse_arms(gff3: Path) -> Iterator[Arm]:
    """Yield every LTR arm in the flanking-LTR track."""
    with gff3.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != ARM_FEATURE:
                continue
            attributes = _attributes(fields[8])
            parent = attributes.get("Parent", "")
            if not parent:
                continue
            yield Arm(
                seqname=fields[0],
                start=int(fields[3]),
                end=int(fields[4]),
                parent=parent,
                arm=attributes.get("arm", "?"),
            )


def select_bait(arms: Iterator[Arm], parents: set[str], min_length: int) -> list[Arm]:
    """Keep arms of ERV-bearing elements that are long enough to be informative."""
    return [arm for arm in arms if arm.parent in parents and arm.length >= min_length]


def write_bed(arms: list[Arm], path: Path) -> int:
    """Write BED6 and return the row count.

    GFF3 is 1-based with both ends inclusive; BED is 0-based half-open, so the start
    moves back one and the end stays. Strand is always `+`: blastn searches both
    strands regardless, and extracting the arm as written keeps its coordinates
    directly comparable to the track it came from.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for arm in arms:
            handle.write(
                f"{arm.seqname}\t{arm.start - 1}\t{arm.end}\t{arm.name}\t.\t+\n"
            )
    return len(arms)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--flanking-ltr-gff3", type=Path, required=True)
    parser.add_argument("--loci-csv", type=Path, required=True)
    parser.add_argument("--out-bed", type=Path, required=True)
    parser.add_argument(
        "--min-bait-length",
        type=int,
        required=True,
        help="From config solo_ltr.min_bait_length. Required on purpose: a default "
        "here would be a second source of truth that drifts from config.yaml.",
    )
    parser.add_argument("--log", type=Path, help="job log (the Snakemake log: path)")
    args = parser.parse_args(argv)
    job_logging(args.log, "solo_bait_builder")

    parents = erv_bearing_parents(args.loci_csv)
    arms = select_bait(
        parse_arms(args.flanking_ltr_gff3), parents, args.min_bait_length
    )
    written = write_bed(arms, args.out_bed)
    if not written:
        logger.warning(
            "no bait: none of %d ERV-bearing elements has an LTR arm of at least %d "
            "bp, so no solo LTR can be found in this genome. Lower "
            "solo_ltr.min_bait_length if it should have some",
            len(parents),
            args.min_bait_length,
        )
        return
    logger.log(
        OK,
        "%s bait arms from %s of %s ERV-bearing elements",
        f"{written:,}",
        f"{len({a.parent for a in arms}):,}",
        f"{len(parents):,}",
    )


if __name__ == "__main__":
    run_main(main)
