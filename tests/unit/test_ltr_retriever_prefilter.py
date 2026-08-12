"""Unit tests for ``workflow/scripts/solo_ltr/ltr_retriever_prefilter.py``.

Covers the pure helpers (``_parse_valid_ranges``, ``_intervals_overlap``,
``_any_overlap``) and the public ``prefilter_scn`` function, which writes both
the retroviral-restricted SCN (Coupling A, ADR-003) and the unfiltered
passthrough from a single read pass.

Two contracts changed when the suffix-array dependency was removed:

* the ``seq-nr -> chromosome`` mapping now comes from the LTRharvest GFF3's
  ``seq_number=`` attribute instead of the index's ``.des`` file, which is zero
  bytes for Antrozous;
* SCN and GFF3 coordinates are both treated as **1-based closed**. The earlier
  code shifted GFF3 starts by ``- 1`` on the belief that the SCN was 0-based;
  reconstructing a real SCN from its GFF3 byte-for-byte proved the two frames
  are identical, so that shift was a systematic 1 bp widening of every valid
  interval. ``test_prefilter_treats_adjacent_ranges_as_non_overlapping`` is the
  discriminating case.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from ltr_retriever_prefilter import (
    _any_overlap,
    _intervals_overlap,
    _parse_valid_ranges,
    prefilter_scn,
)
from ltr_retriever_prefilter import main as prefilter_main


# ---------------------------------------------------------------------
# _parse_valid_ranges
# ---------------------------------------------------------------------
def test_parse_valid_ranges_keeps_one_based_closed_coordinates(
    tmp_path: Path,
) -> None:
    """GFF3 coordinates pass through unchanged, matching the SCN's own frame."""
    gff = tmp_path / "valid.gff3"
    gff.write_text(
        "##gff-version 3\nchr1\tRetroSeek\tERV\t100\t200\t.\t+\t.\tID=erv1\n"
    )
    assert _parse_valid_ranges(gff) == {"chr1": [(100, 200)]}


def test_parse_valid_ranges_skips_comments_and_short_rows(tmp_path: Path) -> None:
    """Comment lines and rows with fewer than 5 fields are ignored silently."""
    gff = tmp_path / "messy.gff3"
    gff.write_text(
        "# leading comment\n"
        "\n"  # blank
        "chr1\ttoo\tshort\n"  # only 3 fields
        "chr1\tRetroSeek\tERV\t10\t20\t.\t+\t.\tID=ok\n"
        "## another comment\n"
    )
    assert _parse_valid_ranges(gff) == {"chr1": [(10, 20)]}


def test_parse_valid_ranges_sorts_by_start(tmp_path: Path) -> None:
    """Per-chromosome intervals are sorted to enable short-circuit overlap walks."""
    gff = tmp_path / "unsorted.gff3"
    gff.write_text(
        "chr1\t.\t.\t500\t600\t.\t.\t.\tID=a\n"
        "chr1\t.\t.\t100\t200\t.\t.\t.\tID=b\n"
        "chr1\t.\t.\t300\t400\t.\t.\t.\tID=c\n"
    )
    assert _parse_valid_ranges(gff)["chr1"] == [(100, 200), (300, 400), (500, 600)]


def test_parse_valid_ranges_missing_file_raises(tmp_path: Path) -> None:
    missing = tmp_path / "no_such.gff3"
    with pytest.raises(FileNotFoundError, match="valid_ranges GFF3 not found"):
        _parse_valid_ranges(missing)


# ---------------------------------------------------------------------
# _intervals_overlap
# ---------------------------------------------------------------------
def test_intervals_overlap_closed_inclusive_boundary() -> None:
    """Touching at a single point counts as overlap (closed-closed semantics).

    ``[5, 10]`` and ``[10, 15]`` share the point 10 -> True. Off by one
    (``[5, 9]`` and ``[10, 15]``) -> False.
    """
    assert _intervals_overlap(5, 10, 10, 15) is True
    assert _intervals_overlap(5, 9, 10, 15) is False


def test_intervals_overlap_full_containment() -> None:
    """An interval contained inside another overlaps."""
    assert _intervals_overlap(50, 60, 0, 100) is True


def test_intervals_overlap_disjoint() -> None:
    assert _intervals_overlap(0, 5, 10, 20) is False


# ---------------------------------------------------------------------
# _any_overlap
# ---------------------------------------------------------------------
def test_any_overlap_returns_false_on_empty_or_none() -> None:
    assert _any_overlap(0, 100, None) is False
    assert _any_overlap(0, 100, []) is False


def test_any_overlap_short_circuits_when_start_exceeds_end() -> None:
    """The walk must abort once an interval's start is past our end.

    Given a sorted list, any later interval's start is also past our
    end - so further iteration is wasted.
    """
    intervals = [(0, 5), (10, 20), (50, 60), (100, 200)]
    assert _any_overlap(7, 9, intervals) is False
    assert _any_overlap(15, 16, intervals) is True


# ---------------------------------------------------------------------
# prefilter_scn - dual-output contract
#
# The prefilter writes both the retroviral-restricted SCN (Coupling A,
# default consumer for LTR_retriever) AND the unfiltered ``_full.scn``
# from a single read pass. Which one feeds LTR_retriever is a runtime
# decision driven by ``config.ltr_retriever.source_scn``, but both are
# always materialised so a user can inspect either at will.
# ---------------------------------------------------------------------
def _write_tiny_fixture(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Build a tiny SCN + LTRharvest GFF3 + valid_ranges trio for prefilter tests.

    The LTRharvest GFF3 supplies only the ``seq-nr -> chromosome`` mapping here
    (seq 0 = chr1, seq 1 = chr2); its element coordinates are irrelevant to the
    prefilter, which reads spans from the SCN itself.
    """
    scn = tmp_path / "tiny.scn"
    scn.write_text(
        "# LTR_FINDER args=...\n"
        "# predictions are reported in the following way\n"
        "# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) "
        "s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr\n"
        "100 500 401 100 200 101 400 500 101 95.0 0\n"
        "1500 1900 401 1500 1600 101 1800 1900 101 92.0 0\n"
        "100 500 401 100 200 101 400 500 101 95.0 1\n"
    )
    harvest = tmp_path / "tiny_ltrharvest.gff3"
    harvest.write_text(
        "chr1\tLTRharvest\tLTR_retrotransposon\t100\t500\t.\t?\t.\t"
        "ID=LTR_retrotransposon1;ltr_similarity=95.00;seq_number=0\n"
        "chr1\tLTRharvest\tlong_terminal_repeat\t100\t200\t.\t?\t.\t"
        "Parent=LTR_retrotransposon1\n"
        "chr1\tLTRharvest\tlong_terminal_repeat\t400\t500\t.\t?\t.\t"
        "Parent=LTR_retrotransposon1\n"
        "chr2\tLTRharvest\tLTR_retrotransposon\t100\t500\t.\t?\t.\t"
        "ID=LTR_retrotransposon2;ltr_similarity=95.00;seq_number=1\n"
        "chr2\tLTRharvest\tlong_terminal_repeat\t100\t200\t.\t?\t.\t"
        "Parent=LTR_retrotransposon2\n"
        "chr2\tLTRharvest\tlong_terminal_repeat\t400\t500\t.\t?\t.\t"
        "Parent=LTR_retrotransposon2\n"
    )
    gff = tmp_path / "tiny_valid.gff3"
    gff.write_text(
        "chr1\tRetroSeek\tERV\t101\t501\t.\t+\t.\tID=erv1\n"
        # No row for chr2 -> seq-nr=1 candidate is in full but not retroviral.
    )
    return scn, harvest, gff


def test_prefilter_writes_both_retroviral_and_full_outputs(tmp_path: Path) -> None:
    """Single invocation produces both files; retroviral is a subset of full."""
    scn, harvest, gff = _write_tiny_fixture(tmp_path)
    retroviral = tmp_path / "tiny_retroviral.scn"
    full = tmp_path / "tiny_full.scn"

    prefilter_scn(scn, harvest, gff, retroviral, full)

    assert retroviral.is_file()
    assert full.is_file()
    retroviral_data = [
        line
        for line in retroviral.read_text().splitlines()
        if not line.startswith("#") and line.strip()
    ]
    full_data = [
        line
        for line in full.read_text().splitlines()
        if not line.startswith("#") and line.strip()
    ]
    assert len(retroviral_data) == 1
    assert len(full_data) == 3
    assert set(retroviral_data).issubset(set(full_data))


def test_prefilter_treats_adjacent_ranges_as_non_overlapping(tmp_path: Path) -> None:
    """The regression for the retired ``start - 1`` shift.

    An SCN candidate spanning [100, 500] and a valid range starting at 501 are
    adjacent, not overlapping. The old code shifted the GFF3 start to 500 and so
    counted them as overlapping, admitting candidates that share no base with
    any domain-validated ERV.
    """
    scn = tmp_path / "adjacent.scn"
    scn.write_text("100 500 401 100 200 101 400 500 101 95.0 0\n")
    harvest = tmp_path / "adjacent_ltrharvest.gff3"
    harvest.write_text(
        "chr1\tLTRharvest\tLTR_retrotransposon\t100\t500\t.\t?\t.\t"
        "ID=LTR_retrotransposon1;ltr_similarity=95.00;seq_number=0\n"
        "chr1\tLTRharvest\tlong_terminal_repeat\t100\t200\t.\t?\t.\t"
        "Parent=LTR_retrotransposon1\n"
        "chr1\tLTRharvest\tlong_terminal_repeat\t400\t500\t.\t?\t.\t"
        "Parent=LTR_retrotransposon1\n"
    )
    gff = tmp_path / "adjacent_valid.gff3"
    gff.write_text("chr1\tRetroSeek\tERV\t501\t600\t.\t+\t.\tID=erv1\n")

    rows_in, kept_retroviral, kept_full = prefilter_scn(
        scn, harvest, gff, tmp_path / "r.scn", tmp_path / "f.scn"
    )
    assert (rows_in, kept_retroviral, kept_full) == (1, 0, 1)


def test_prefilter_keeps_a_candidate_sharing_a_single_base(tmp_path: Path) -> None:
    """The other side of the boundary: one shared base is a real overlap."""
    scn = tmp_path / "touching.scn"
    scn.write_text("100 500 401 100 200 101 400 500 101 95.0 0\n")
    harvest = tmp_path / "touching_ltrharvest.gff3"
    harvest.write_text(
        "chr1\tLTRharvest\tLTR_retrotransposon\t100\t500\t.\t?\t.\t"
        "ID=LTR_retrotransposon1;ltr_similarity=95.00;seq_number=0\n"
        "chr1\tLTRharvest\tlong_terminal_repeat\t100\t200\t.\t?\t.\t"
        "Parent=LTR_retrotransposon1\n"
        "chr1\tLTRharvest\tlong_terminal_repeat\t400\t500\t.\t?\t.\t"
        "Parent=LTR_retrotransposon1\n"
    )
    gff = tmp_path / "touching_valid.gff3"
    gff.write_text("chr1\tRetroSeek\tERV\t500\t600\t.\t+\t.\tID=erv1\n")

    _, kept_retroviral, _ = prefilter_scn(
        scn, harvest, gff, tmp_path / "r.scn", tmp_path / "f.scn"
    )
    assert kept_retroviral == 1


def test_prefilter_full_output_byte_equivalent_to_input(tmp_path: Path) -> None:
    """Without malformed rows, ``_full.scn`` is byte-equal to the source SCN.

    This is the load-bearing invariant for ``source_scn: full`` mode:
    LTR_retriever fed the full SCN must see exactly what LTRharvest
    emitted, including comments. Locks against any future regression
    where the prefilter accidentally rewrites or reformats data rows.
    """
    scn, harvest, gff = _write_tiny_fixture(tmp_path)
    retroviral = tmp_path / "out_retroviral.scn"
    full = tmp_path / "out_full.scn"

    prefilter_scn(scn, harvest, gff, retroviral, full)

    assert full.read_bytes() == scn.read_bytes()


def test_prefilter_preserves_header_comments_in_both_outputs(tmp_path: Path) -> None:
    """Both SCN files preserve LTRharvest's header comments verbatim."""
    scn, harvest, gff = _write_tiny_fixture(tmp_path)
    retroviral = tmp_path / "r.scn"
    full = tmp_path / "f.scn"
    prefilter_scn(scn, harvest, gff, retroviral, full)
    for output in (retroviral, full):
        text = output.read_text()
        assert "# LTR_FINDER args=..." in text
        assert "# s(ret) e(ret)" in text


def test_prefilter_returns_three_counts_tuple(tmp_path: Path) -> None:
    """Return shape: ``(rows_in, rows_kept_retroviral, rows_kept_full)``.

    ``rows_kept_full`` always equals ``rows_in`` modulo malformed rows
    (both in this test are zero, so equality holds).
    """
    scn, harvest, gff = _write_tiny_fixture(tmp_path)
    result = prefilter_scn(scn, harvest, gff, tmp_path / "r.scn", tmp_path / "f.scn")
    assert result == (3, 1, 3)


def test_prefilter_drops_rows_whose_seq_nr_is_unknown(tmp_path: Path) -> None:
    """A seq-nr absent from the GFF3 cannot be mapped, so it is not retroviral.

    It still reaches the full passthrough, which must mirror LTRharvest exactly.
    """
    scn = tmp_path / "unknown.scn"
    scn.write_text("100 500 401 100 200 101 400 500 101 95.0 7\n")
    harvest = tmp_path / "unknown_ltrharvest.gff3"
    harvest.write_text(
        "chr1\tLTRharvest\tLTR_retrotransposon\t100\t500\t.\t?\t.\t"
        "ID=LTR_retrotransposon1;ltr_similarity=95.00;seq_number=0\n"
        "chr1\tLTRharvest\tlong_terminal_repeat\t100\t200\t.\t?\t.\t"
        "Parent=LTR_retrotransposon1\n"
        "chr1\tLTRharvest\tlong_terminal_repeat\t400\t500\t.\t?\t.\t"
        "Parent=LTR_retrotransposon1\n"
    )
    gff = tmp_path / "unknown_valid.gff3"
    gff.write_text("chr1\tRetroSeek\tERV\t100\t600\t.\t+\t.\tID=erv1\n")

    rows_in, kept_retroviral, kept_full = prefilter_scn(
        scn, harvest, gff, tmp_path / "r.scn", tmp_path / "f.scn"
    )
    assert (rows_in, kept_retroviral, kept_full) == (1, 0, 1)


def test_prefilter_main_cli_accepts_dual_output_flags(tmp_path: Path) -> None:
    """``--output-retroviral`` and ``--output-full`` are both required CLI args."""
    scn, harvest, gff = _write_tiny_fixture(tmp_path)
    retroviral = tmp_path / "cli_retroviral.scn"
    full = tmp_path / "cli_full.scn"

    rc = prefilter_main(
        [
            "--scn",
            str(scn),
            "--ltrharvest-gff3",
            str(harvest),
            "--valid-ranges",
            str(gff),
            "--output-retroviral",
            str(retroviral),
            "--output-full",
            str(full),
        ]
    )
    assert rc == 0
    assert retroviral.is_file()
    assert full.is_file()
