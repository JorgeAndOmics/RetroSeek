"""Pre-filter LTRharvest SCN output by intersecting with valid_ranges.gff3.

This is *Coupling A* of the LTR_retriever integration. RetroSeek's
``ranges_analysis`` step has already validated a subset of LTRharvest
candidates as retroviral by domain matching; this script writes two SCN
files from a single read pass:

* ``{genome}_retroviral.scn`` - the rows whose paired-LTR coordinates
  overlap a valid_ranges interval on the same chromosome. This is the
  default LTR_retriever input under ``config.ltr_retriever.source_scn:
  retroviral``: LTR_retriever's family-building and BLAST-back passes
  see only retroviral-confirmed candidates, so the solo LTRs it
  discovers are guaranteed retroviral.
* ``{genome}_full.scn`` - every well-formed row from the source SCN,
  comments included, byte-equivalent to the source modulo malformed
  rows. This is the LTR_retriever input under
  ``config.ltr_retriever.source_scn: full``, useful for non-retroviral
  exploration or cross-validation against the retroviral output.

Both files always materialise; the runtime decision lives one rule
downstream in ``run_ltr_retriever.py``.

Input files
-----------
LTRharvest SCN
    The screen-output format LTRharvest writes to stdout. One row per
    detected paired-LTR structure with columns:

        s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr

    The ``seq-nr`` column is the LTRharvest-internal sequence index (0,
    1, 2, ...) into the input FASTA, NOT a chromosome name.

LTRharvest GFF3
    The same predictions in GFF3 form. Read here only for its
    ``seq_number=`` attribute, which sits on the same row as the
    chromosome name and so supplies the ``seq-nr -> chromosome``
    mapping. This replaces the suffix array's ``.des`` file: Antrozous's
    index files are all zero bytes, and rebuilding the index would
    cascade LTRharvest and LTRdigest across every genome (see
    ``scn_from_ltrharvest_gff3.py`` for the full reasoning).

valid_ranges.gff3
    RetroSeek's domain-validated retroviral ERV track - output of
    ``ranges_analysis_setup``. Standard GFF3: comment lines start with
    ``#``, data rows have 9 tab-separated fields.

Coordinate systems
------------------
Both the SCN and GFF3 are **1-based closed**, so coordinates compare
directly with no shift. This was verified empirically by reconstructing a
real SCN from its GFF3: ``s(ret)`` equals the GFF3 start exactly (see the
byte-equality regression in ``tests/unit/test_scn_from_ltrharvest_gff3.py``).
An earlier revision believed the SCN was 0-based and normalised GFF3 starts
via ``start - 1``, which widened every valid interval by one base and
admitted candidates sharing no base with any validated ERV.

Two intervals overlap iff ``a_start <= b_end`` AND ``a_end >= b_start``.

Usage (CLI)
-----------
::

    python ltr_retriever_prefilter.py \\
        --scn data/ltr_scn/{genome}.scn \\
        --ltrharvest-gff3 results/tracks/ltrharvest/{genome}.gff3 \\
        --valid-ranges results/tracks/valid/{genome}.gff3 \\
        --output-retroviral data/ltr_scn/{genome}_retroviral.scn \\
        --output-full       data/ltr_scn/{genome}_full.scn
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections import defaultdict
from pathlib import Path

from scn_from_ltrharvest_gff3 import parse_ltrharvest_gff3, seq_number_map

logger = logging.getLogger(__name__)


def _parse_valid_ranges(gff3_path: Path) -> dict[str, list[tuple[int, int]]]:
    """Return per-chromosome sorted list of (start, end) intervals, 1-based closed.

    GFF3 coordinates pass through unchanged: they share the SCN's frame, so
    no conversion is needed for the overlap test.

    Comments (``#``-prefixed lines) and malformed rows (< 9 fields) are
    silently skipped.
    """
    if not gff3_path.exists():
        raise FileNotFoundError(f"valid_ranges GFF3 not found: {gff3_path}")
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with gff3_path.open() as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            try:
                seqid = fields[0]
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                continue
            intervals[seqid].append((start, end))
    # Sort per-chromosome so overlap checks can short-circuit.
    for chrom_intervals in intervals.values():
        chrom_intervals.sort()
    return dict(intervals)


def _intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """Closed-interval overlap test. True iff the intersection is non-empty."""
    return a_start <= b_end and a_end >= b_start


def _any_overlap(start: int, end: int, intervals: list[tuple[int, int]] | None) -> bool:
    """Return True if ``[start, end]`` overlaps any interval in the sorted list."""
    if not intervals:
        return False
    # Sorted by start; walk until a start exceeds our end.
    for int_start, int_end in intervals:
        if int_start > end:
            return False
        if _intervals_overlap(start, end, int_start, int_end):
            return True
    return False


def prefilter_scn(
    scn_path: Path,
    ltrharvest_gff3_path: Path,
    valid_ranges_path: Path,
    retroviral_output_path: Path,
    full_output_path: Path,
) -> tuple[int, int, int]:
    """Filter SCN rows in a single pass; write retroviral + full SCN files.

    Comments (``#``-prefixed) and blank lines are written verbatim to
    both outputs so LTR_retriever's parser sees a well-formed SCN
    regardless of which path is selected downstream.

    Each well-formed data row is always written to ``full_output_path``;
    it is *also* written to ``retroviral_output_path`` iff its paired-LTR
    coordinates overlap a valid_ranges interval on the same chromosome.

    Returns
    -------
    tuple[int, int, int]
        ``(rows_in, rows_kept_retroviral, rows_kept_full)`` - input row
        count, the count retained by the retroviral filter, and the
        count emitted to the full output. ``rows_kept_full`` equals
        ``rows_in`` modulo malformed rows.
    """
    chrom_names = seq_number_map(parse_ltrharvest_gff3(ltrharvest_gff3_path))
    valid_intervals = _parse_valid_ranges(valid_ranges_path)

    rows_in = 0
    rows_kept_retroviral = 0
    rows_kept_full = 0

    retroviral_output_path.parent.mkdir(parents=True, exist_ok=True)
    full_output_path.parent.mkdir(parents=True, exist_ok=True)

    with (
        scn_path.open() as fin,
        retroviral_output_path.open("w") as fout_retroviral,
        full_output_path.open("w") as fout_full,
    ):
        for raw in fin:
            # Preserve SCN header comment lines verbatim - LTR_retriever parses them.
            if raw.startswith("#") or not raw.strip():
                fout_retroviral.write(raw)
                fout_full.write(raw)
                continue
            parts = raw.split()
            if len(parts) < 11:
                # Malformed row - skip (LTR_retriever would likely skip it too).
                continue
            rows_in += 1
            try:
                ret_start = int(parts[0])
                ret_end = int(parts[1])
                seq_nr = int(parts[10])
            except ValueError:
                continue
            # Always emit to the full file - that's the load-bearing
            # contract for ``source_scn: full`` mode.
            fout_full.write(raw)
            rows_kept_full += 1
            chrom = chrom_names.get(seq_nr)
            if chrom is None:
                # seq_nr references a sequence the GFF3 never mentions - drop
                # from retroviral, but it is already in full.
                continue
            if _any_overlap(ret_start, ret_end, valid_intervals.get(chrom)):
                fout_retroviral.write(raw)
                rows_kept_retroviral += 1

    return rows_in, rows_kept_retroviral, rows_kept_full


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--scn", type=Path, required=True, help="Input LTRharvest SCN file."
    )
    parser.add_argument(
        "--ltrharvest-gff3",
        type=Path,
        required=True,
        help="LTRharvest GFF3; supplies the seq-nr -> chromosome mapping.",
    )
    parser.add_argument(
        "--valid-ranges",
        type=Path,
        required=True,
        help="Path to valid_ranges.gff3 from ranges_analysis_setup.",
    )
    parser.add_argument(
        "--output-retroviral",
        type=Path,
        required=True,
        help="Output path for the retroviral-restricted SCN (Coupling A).",
    )
    parser.add_argument(
        "--output-full",
        type=Path,
        required=True,
        help="Output path for the full unfiltered SCN (passthrough copy).",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    logger.info("Pre-filtering SCN %s against %s", args.scn, args.valid_ranges)

    rows_in, rows_kept_retroviral, rows_kept_full = prefilter_scn(
        scn_path=args.scn,
        ltrharvest_gff3_path=args.ltrharvest_gff3,
        valid_ranges_path=args.valid_ranges,
        retroviral_output_path=args.output_retroviral,
        full_output_path=args.output_full,
    )
    pct_retroviral = (100.0 * rows_kept_retroviral / rows_in) if rows_in else 0.0
    logger.info(
        "Retained %d of %d SCN rows (%.1f%%) for retroviral; "
        "wrote %d rows to full output",
        rows_kept_retroviral,
        rows_in,
        pct_retroviral,
        rows_kept_full,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
