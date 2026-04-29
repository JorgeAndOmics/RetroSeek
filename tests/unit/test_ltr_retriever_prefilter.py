"""Unit tests for ``workflow/scripts/ltr_retriever_prefilter.py``.

Covers the four pure helpers (``_parse_des``, ``_parse_valid_ranges``,
``_intervals_overlap``, ``_any_overlap``) and the public ``prefilter_scn``
function. Tests pin the current single-output contract; Phase 2 extends
these with dual-output (retroviral + full SCN) cases.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from ltr_retriever_prefilter import (
    _any_overlap,
    _intervals_overlap,
    _parse_des,
    _parse_valid_ranges,
    prefilter_scn,
)


# ---------------------------------------------------------------------
# _parse_des
# ---------------------------------------------------------------------
def test_parse_des_basic_text_file(tmp_path: Path) -> None:
    """Plain text .des file: one chromosome name per line, returned in order."""
    des = tmp_path / "basic.des"
    des.write_bytes(b"chr1\nchr2\nchr3\n")
    assert _parse_des(des) == ["chr1", "chr2", "chr3"]


def test_parse_des_handles_binary_trailer(tmp_path: Path) -> None:
    """Real-world .des files end with a 0xFF-padded trailer.

    Reproduces the failure mode that motivated the binary-mode reader:
    a UTF-8 text reader chokes on 0xFF mid-stream. The parser must
    treat the first 0xFF byte as end-of-data and return only the
    chromosome names that preceded it.
    """
    des = tmp_path / "with_trailer.des"
    des.write_bytes(b"chr1\nchr2\nchr3\n" + b"\xff" * 16)
    assert _parse_des(des) == ["chr1", "chr2", "chr3"]


def test_parse_des_strips_to_first_token(tmp_path: Path) -> None:
    """LTRharvest -seqids uses only the first whitespace-delimited token.

    A full FASTA header like ``chr1 dna:chromosome chromosome:GRCh38:1``
    must be reduced to ``chr1`` so the chromosome lookup matches the
    GFF3 seqid column.
    """
    des = tmp_path / "with_descriptions.des"
    des.write_bytes(b"chr1 dna:chromosome chromosome:GRCh38:1\nchr2 some other text\n")
    assert _parse_des(des) == ["chr1", "chr2"]


def test_parse_des_missing_file_raises(tmp_path: Path) -> None:
    """An absent .des file must raise FileNotFoundError, not silently empty out."""
    missing = tmp_path / "does_not_exist.des"
    with pytest.raises(FileNotFoundError, match="descriptor file not found"):
        _parse_des(missing)


# ---------------------------------------------------------------------
# _parse_valid_ranges
# ---------------------------------------------------------------------
def test_parse_valid_ranges_converts_one_indexed_to_zero_indexed(tmp_path: Path) -> None:
    """GFF3 1-indexed closed → returned as 0-indexed (start - 1, end unchanged).

    Single row at GFF3 [100, 200] inclusive should become [99, 200] in
    the 0-indexed half-open / closed convention used by the SCN overlap
    checks.
    """
    gff = tmp_path / "valid.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tRetroSeek\tERV\t100\t200\t.\t+\t.\tID=erv1\n"
    )
    intervals = _parse_valid_ranges(gff)
    assert intervals == {"chr1": [(99, 200)]}


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
    intervals = _parse_valid_ranges(gff)
    assert intervals == {"chr1": [(9, 20)]}


def test_parse_valid_ranges_sorts_by_start(tmp_path: Path) -> None:
    """Per-chromosome intervals are sorted to enable short-circuit overlap walks."""
    gff = tmp_path / "unsorted.gff3"
    gff.write_text(
        "chr1\t.\t.\t500\t600\t.\t.\t.\tID=a\n"
        "chr1\t.\t.\t100\t200\t.\t.\t.\tID=b\n"
        "chr1\t.\t.\t300\t400\t.\t.\t.\tID=c\n"
    )
    intervals = _parse_valid_ranges(gff)
    assert intervals["chr1"] == [(99, 200), (299, 400), (499, 600)]


def test_parse_valid_ranges_missing_file_raises(tmp_path: Path) -> None:
    missing = tmp_path / "no_such.gff3"
    with pytest.raises(FileNotFoundError, match="valid_ranges GFF3 not found"):
        _parse_valid_ranges(missing)


# ---------------------------------------------------------------------
# _intervals_overlap
# ---------------------------------------------------------------------
def test_intervals_overlap_closed_inclusive_boundary() -> None:
    """Touching at a single point counts as overlap (closed-closed semantics).

    ``[5, 10]`` and ``[10, 15]`` share the point 10 → True. Off by one
    (``[5, 9]`` and ``[10, 15]``) → False.
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
    end — so further iteration is wasted.
    """
    intervals = [(0, 5), (10, 20), (50, 60), (100, 200)]
    assert _any_overlap(7, 9, intervals) is False
    assert _any_overlap(15, 16, intervals) is True


# ---------------------------------------------------------------------
# prefilter_scn (current single-output contract; Phase 2 will extend)
# ---------------------------------------------------------------------
def _write_tiny_fixture(tmp_path: Path) -> tuple[Path, Path, Path]:
    """Build a tiny SCN + .des + valid_ranges trio for prefilter tests."""
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
    des = tmp_path / "tiny.des"
    des.write_bytes(b"chr1\nchr2\n")
    gff = tmp_path / "tiny_valid.gff3"
    gff.write_text(
        "chr1\tRetroSeek\tERV\t101\t501\t.\t+\t.\tID=erv1\n"
        # No row for chr2 → seq-nr=1 candidate must be filtered out.
    )
    return scn, des, gff


def test_prefilter_scn_keeps_only_overlapping_rows(tmp_path: Path) -> None:
    scn, des, gff = _write_tiny_fixture(tmp_path)
    out = tmp_path / "out.scn"
    rows_in, rows_kept = prefilter_scn(scn, des, gff, out)
    assert rows_in == 3
    assert rows_kept == 1


def test_prefilter_scn_preserves_header_comments(tmp_path: Path) -> None:
    scn, des, gff = _write_tiny_fixture(tmp_path)
    out = tmp_path / "out.scn"
    prefilter_scn(scn, des, gff, out)
    output_text = out.read_text()
    assert "# LTR_FINDER args=..." in output_text
    assert "# s(ret) e(ret)" in output_text


def test_prefilter_scn_returns_correct_counts(tmp_path: Path) -> None:
    """Return tuple matches the (rows_in, rows_kept) contract."""
    scn, des, gff = _write_tiny_fixture(tmp_path)
    out = tmp_path / "out.scn"
    result = prefilter_scn(scn, des, gff, out)
    assert isinstance(result, tuple)
    assert len(result) == 2
    assert result == (3, 1)
