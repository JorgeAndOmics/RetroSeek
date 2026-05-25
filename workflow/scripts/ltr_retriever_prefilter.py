"""Pre-filter LTRharvest SCN output by intersecting with valid_ranges.gff3.

This is *Coupling A* of the LTR_retriever integration. RetroSeek's
``ranges_analysis`` step has already validated a subset of LTRharvest
candidates as retroviral by domain matching; this script writes two SCN
files from a single read pass:

* ``{genome}_retroviral.scn`` — the rows whose paired-LTR coordinates
  overlap a valid_ranges interval on the same chromosome. This is the
  default LTR_retriever input under ``config.ltr_retriever.source_scn:
  retroviral``: LTR_retriever's family-building and BLAST-back passes
  see only retroviral-confirmed candidates, so the solo LTRs it
  discovers are guaranteed retroviral.
* ``{genome}_full.scn`` — every well-formed row from the source SCN,
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
    detected paired-LTR structure with columns (1-indexed per the tool's
    convention, but positions are in the 0-indexed subject sequence
    frame used by LTRharvest):

        s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr

    The ``seq-nr`` column is the LTRharvest-internal sequence index (0,
    1, 2, ...) into the input FASTA, NOT a chromosome name.

``.des`` file
    Companion to the suffix-array index. Line ``N`` (0-indexed) holds
    the chromosome name assigned to ``seq-nr = N``. Used to translate
    the SCN's ``seq-nr`` column into chromosome names that match
    valid_ranges.gff3.

valid_ranges.gff3
    RetroSeek's domain-validated retroviral ERV track — output of
    ``ranges_analysis_setup``. Standard GFF3: comment lines start with
    ``#``, data rows have 9 tab-separated fields.

Coordinate systems
------------------
LTRharvest SCN ``s(ret)`` / ``e(ret)`` are 0-indexed closed intervals.
GFF3 uses 1-indexed closed intervals. The overlap check normalises
GFF3 to 0-indexed closed via ``start - 1``; the SCN values are left
as-is. Two intervals overlap iff ``a_start <= b_end`` AND
``a_end >= b_start``.

Usage (CLI)
-----------
::

    python ltr_retriever_prefilter.py \\
        --scn data/ltr_scn/{genome}.scn \\
        --des SPECIES_DB/{genome}/{genome}.des \\
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

logger = logging.getLogger(__name__)


def _parse_des(des_path: Path) -> list[str]:
    """Return the list of chromosome names indexed by LTRharvest seq-nr.

    The ``.des`` file produced by ``gt suffixerator`` contains one
    FASTA-header-like string per line, in the order LTRharvest numbers
    them (``seq-nr=0`` → first line). Whitespace is trimmed; empty
    lines are skipped. Only the first whitespace-separated token is
    retained since LTRharvest uses that as the seqid when ``-seqids``
    is active.

    Raises
    ------
    FileNotFoundError
        If ``des_path`` does not exist.
    """
    if not des_path.exists():
        raise FileNotFoundError(f"LTRharvest descriptor file not found: {des_path}")
    # gt suffixerator's `.des` is mostly text (one FASTA header per `\n`-
    # delimited line) but ends with a binary 0xFF-padded trailer. Reading
    # it as text crashes UTF-8 decoding mid-stream. Open binary, split on
    # `\n`, decode each line leniently, and stop when we hit non-printable
    # content (the trailer's 0xFF bytes are a reliable end-of-data sentinel).
    names: list[str] = []
    with des_path.open("rb") as handle:
        raw_bytes = handle.read()
    for raw_line in raw_bytes.split(b"\n"):
        if not raw_line:
            continue
        # Trailer detection: a line that is all 0xFF (or starts with 0xFF) is
        # the binary padding, not a chromosome name.
        if raw_line[0] == 0xFF:
            break
        try:
            line = raw_line.decode("utf-8").strip()
        except UnicodeDecodeError:
            # Mixed-codec garbage — treat as end of textual content.
            break
        if not line:
            continue
        # Keep only the first token — matches LTRharvest -seqids behaviour.
        names.append(line.split()[0])
    return names


def _parse_valid_ranges(gff3_path: Path) -> dict[str, list[tuple[int, int]]]:
    """Return per-chromosome sorted list of (start, end) intervals, 0-indexed.

    GFF3 is 1-indexed closed; this converts start to 0-indexed (``- 1``)
    and leaves end as-is (0-indexed closed ↔ 1-indexed closed have the
    same right edge under this convention).

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
                start = int(fields[3]) - 1  # GFF3 1-indexed → 0-indexed
                end = int(fields[4])
            except ValueError:
                continue
            intervals[seqid].append((start, end))
    # Sort per-chromosome so overlap checks can short-circuit.
    for chrom_intervals in intervals.values():
        chrom_intervals.sort()
    return dict(intervals)


def _intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """Closed-interval overlap test (0-indexed). True iff a ∩ b is non-empty."""
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
    des_path: Path,
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
        ``(rows_in, rows_kept_retroviral, rows_kept_full)`` — input row
        count, the count retained by the retroviral filter, and the
        count emitted to the full output. ``rows_kept_full`` equals
        ``rows_in`` modulo malformed rows.
    """
    chrom_names = _parse_des(des_path)
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
            # Preserve SCN header comment lines verbatim — LTR_retriever parses them.
            if raw.startswith("#") or not raw.strip():
                fout_retroviral.write(raw)
                fout_full.write(raw)
                continue
            parts = raw.split()
            if len(parts) < 11:
                # Malformed row — skip (LTR_retriever would likely skip it too).
                continue
            rows_in += 1
            try:
                ret_start = int(parts[0])
                ret_end = int(parts[1])
                seq_nr = int(parts[10])
            except ValueError:
                continue
            # Always emit to the full file — that's the load-bearing
            # contract for ``source_scn: full`` mode.
            fout_full.write(raw)
            rows_kept_full += 1
            if not 0 <= seq_nr < len(chrom_names):
                # seq_nr references a chromosome we don't know about — drop
                # from retroviral but it's already in full.
                continue
            chrom = chrom_names[seq_nr]
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
        "--des",
        type=Path,
        required=True,
        help="LTRharvest/suffixerator .des file (one chromosome name per line).",
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
        des_path=args.des,
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
