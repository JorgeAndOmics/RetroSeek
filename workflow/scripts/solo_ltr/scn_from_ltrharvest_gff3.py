"""Reconstruct an LTRharvest SCN file from the GFF3 LTRharvest already wrote.

Why this exists
---------------
``gt ltrharvest`` writes its screen-format predictions (``.scn``) to stdout and
an equivalent GFF3 to a file. The SCN is LTR_retriever's primary input
(``-inharvest``), but it is the more fragile of the two artefacts: it is a
stdout redirect, so a truncated or lost redirect leaves the GFF3 intact and the
SCN gone.

That is exactly the state Antrozous pallidus is in - no ``.scn``, and every
suffix-array file (``.des``, ``.esq``, ``.lcp``, ``.llv``, ``.suf``) zero bytes
after a space-reclaim event. Regenerating the index to re-run LTRharvest would
cost ~20 GB and many hours, and would *also* invalidate LTRharvest and LTRdigest
for every other genome, because ``ltr_harvester_setup`` declares
``input: rules.ltr_index_generator.input`` - an ``expand()`` over all species,
not just the one being harvested. One repair would cascade into a full-pipeline
re-run.

None of that is necessary. The SCN carries no information the GFF3 lacks:

===========================  ==========================================
SCN column                   GFF3 source
===========================  ==========================================
``s(ret)`` / ``e(ret)``      ``LTR_retrotransposon`` start / end
``l(ret)``                   computed, ``end - start + 1``
``s(lLTR)`` / ``e(lLTR)``    first ``long_terminal_repeat`` child
``s(rLTR)`` / ``e(rLTR)``    second ``long_terminal_repeat`` child
``sim(LTRs)``                ``ltr_similarity=`` attribute
``seq-nr``                   ``seq_number=`` attribute
===========================  ==========================================

Verified on Desmodus rotundus, the one model genome holding both artefacts:
all 9,893 reconstructed data rows are byte-identical to the real SCN. See
``tests/unit/test_scn_from_ltrharvest_gff3.py``.

The ``seq_number=`` attribute sits on the same GFF3 row as the chromosome name,
so this module also supplies the ``seq-nr -> chromosome`` mapping that
``ltr_retriever_prefilter.py`` used to read from the suffix array's ``.des``
file. That removes the prefilter's last dependency on the index.

Coordinates
-----------
GFF3 is 1-based closed. The reconstruction copies those values through
unchanged, and the byte-equality regression confirms LTRharvest's SCN uses the
same frame - ``s(ret)`` equals the GFF3 start exactly, not ``start - 1``.

Usage (CLI)
-----------
::

    python scn_from_ltrharvest_gff3.py \\
        --ltrharvest-gff3 results/tracks/ltrharvest/{genome}.gff3 \\
        --output-scn      data/ltr_scn/{genome}.scn
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)

# LTRharvest's own SCN legend, minus the `# args=...` line, which records the
# absolute paths of the original run and cannot be reproduced. LTR_retriever
# skips every `#` line, so this block is documentation for humans reading the
# file directly.
_SCN_HEADER = (
    "# reconstructed by RetroSeek from {source}",
    "# predictions are reported in the following way",
    "# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) "
    "s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr",
    "# where:",
    "# s = starting position",
    "# e = ending position",
    "# l = length",
    "# ret = LTR-retrotransposon",
    "# lLTR = left LTR",
    "# rLTR = right LTR",
    "# sim = similarity",
    "# seq-nr = sequence number",
)

_SIMILARITY_RE = re.compile(r"(?:^|;)ltr_similarity=([0-9.]+)")
_SEQ_NUMBER_RE = re.compile(r"(?:^|;)seq_number=(\d+)")
_ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")


@dataclass(frozen=True)
class LtrharvestRecord:
    """One LTRharvest paired-LTR prediction, in GFF3's 1-based closed frame."""

    seqname: str
    seq_number: int
    ret_start: int
    ret_end: int
    left_start: int
    left_end: int
    right_start: int
    right_end: int
    # Kept as text, not float: the SCN prints two decimals ("90.90"), and
    # round-tripping through float would render that as "90.9" and break
    # byte-equality with LTRharvest's own output.
    similarity: str


def parse_ltrharvest_gff3(gff3_path: Path) -> list[LtrharvestRecord]:
    """Return every paired-LTR prediction in an LTRharvest GFF3, in file order.

    Each ``LTR_retrotransposon`` row opens an element; the two
    ``long_terminal_repeat`` rows that follow are its left and right arms.
    ``repeat_region`` and ``target_site_duplication`` rows are interleaved by
    LTRharvest and ignored here, so arms are matched by which element is open
    rather than by row adjacency.

    Raises
    ------
    FileNotFoundError
        If ``gff3_path`` does not exist.
    ValueError
        If an element does not have exactly two LTR arms, or is missing the
        ``ltr_similarity`` / ``seq_number`` attributes. Both mean a corrupt or
        unexpected GFF3, and an SCN silently short of rows is indistinguishable
        from a genome with fewer LTR elements.
    """
    if not gff3_path.exists():
        raise FileNotFoundError(f"LTRharvest GFF3 not found: {gff3_path}")

    records: list[LtrharvestRecord] = []
    open_id: str | None = None
    open_fields: tuple[str, int, int, int, str] | None = None
    arms: list[tuple[int, int]] = []

    def close_element() -> None:
        """Emit the currently open element, or fail if its arms are wrong."""
        if open_fields is None:
            return
        if len(arms) != 2:
            raise ValueError(
                f"{gff3_path}: element {open_id} has {len(arms)} LTR arm(s), "
                "expected exactly 2. The GFF3 looks truncated or malformed."
            )
        seqname, seq_number, ret_start, ret_end, similarity = open_fields
        (left_start, left_end), (right_start, right_end) = arms
        records.append(
            LtrharvestRecord(
                seqname=seqname,
                seq_number=seq_number,
                ret_start=ret_start,
                ret_end=ret_end,
                left_start=left_start,
                left_end=left_end,
                right_start=right_start,
                right_end=right_end,
                similarity=similarity,
            )
        )

    with gff3_path.open() as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature, attrs = fields[2], fields[8]

            if feature == "LTR_retrotransposon":
                close_element()
                similarity_match = _SIMILARITY_RE.search(attrs)
                seq_number_match = _SEQ_NUMBER_RE.search(attrs)
                id_match = _ID_RE.search(attrs)
                open_id = id_match.group(1) if id_match else "<unnamed>"
                if similarity_match is None or seq_number_match is None:
                    raise ValueError(
                        f"{gff3_path}: element {open_id} is missing "
                        "ltr_similarity= or seq_number=; cannot rebuild the SCN."
                    )
                open_fields = (
                    fields[0],
                    int(seq_number_match.group(1)),
                    int(fields[3]),
                    int(fields[4]),
                    similarity_match.group(1),
                )
                arms = []
            elif feature == "long_terminal_repeat" and open_fields is not None:
                arms.append((int(fields[3]), int(fields[4])))

    close_element()
    return records


def format_scn_row(record: LtrharvestRecord) -> str:
    """Render one record as an LTRharvest SCN data row.

    LTRharvest separates the eleven fields with two spaces and computes each
    length as a closed-interval width (``end - start + 1``). Both are pinned by
    the byte-equality regression, so neither is cosmetic.
    """
    fields = (
        record.ret_start,
        record.ret_end,
        record.ret_end - record.ret_start + 1,
        record.left_start,
        record.left_end,
        record.left_end - record.left_start + 1,
        record.right_start,
        record.right_end,
        record.right_end - record.right_start + 1,
        record.similarity,
        record.seq_number,
    )
    return "  ".join(str(field) for field in fields)


def seq_number_map(records: list[LtrharvestRecord]) -> dict[int, str]:
    """Return the ``seq-nr -> chromosome name`` mapping the SCN needs.

    The SCN's last column is LTRharvest's internal sequence index, not a
    chromosome name. This mapping is what ``ltr_retriever_prefilter.py`` uses to
    compare SCN rows against a GFF3 track, replacing the suffix array's ``.des``
    file (zero bytes for Antrozous).
    """
    return {record.seq_number: record.seqname for record in records}


def write_scn(
    records: list[LtrharvestRecord], output_path: Path, source_gff3: Path
) -> int:
    """Write records as an SCN file (comment header + data rows). Returns the count."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        for line in _SCN_HEADER:
            handle.write(line.format(source=source_gff3) + "\n")
        for record in records:
            handle.write(format_scn_row(record) + "\n")
    return len(records)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--ltrharvest-gff3",
        type=Path,
        required=True,
        help="LTRharvest GFF3 for this genome (results/tracks/ltrharvest/{genome}.gff3).",
    )
    parser.add_argument(
        "--output-scn",
        type=Path,
        required=True,
        help="Output path for the reconstructed SCN (data/ltr_scn/{genome}.scn).",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    records = parse_ltrharvest_gff3(args.ltrharvest_gff3)
    written = write_scn(records, args.output_scn, args.ltrharvest_gff3)
    logger.info(
        "Reconstructed %d SCN rows across %d sequence(s) from %s -> %s",
        written,
        len(seq_number_map(records)),
        args.ltrharvest_gff3,
        args.output_scn,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
