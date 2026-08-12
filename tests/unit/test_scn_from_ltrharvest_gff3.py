"""Unit tests for ``workflow/scripts/solo_ltr/scn_from_ltrharvest_gff3.py``.

The SCN file is a re-rendering of data LTRharvest already wrote to its GFF3:
the retrotransposon span, its two LTR children, ``ltr_similarity=`` and
``seq_number=``. Reconstructing it from the GFF3 lets a genome whose suffix
array is gone (Antrozous, whose index files are all zero bytes) reach
LTR_retriever without regenerating a ~20 GB index - which would also cascade
LTRharvest and LTRdigest across every genome, because ``ltr_harvester_setup``
takes *all* genomes' index files as input.

The load-bearing test here is byte-equality of the data rows against a real
LTRharvest SCN. It is opt-in via ``RETROSEEK_DEV_CACHE`` so no machine path
is committed; the synthetic cases below pin the same contract without it.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest
from scn_from_ltrharvest_gff3 import (
    LtrharvestRecord,
    format_scn_row,
    parse_ltrharvest_gff3,
    seq_number_map,
    write_scn,
)
from scn_from_ltrharvest_gff3 import main as scn_main

# A two-element GFF3 mirroring what `gt ltrharvest -seqids -gff3` emits: each
# LTR_retrotransposon is followed by its two long_terminal_repeat children, and
# the repeat_region / target_site_duplication features are noise for our purpose.
TWO_ELEMENT_GFF3 = """\
##gff-version 3
##sequence-region   chr1 1 1000000
chr1\tLTRharvest\trepeat_region\t24900\t27982\t.\t?\t.\tID=repeat_region1
chr1\tLTRharvest\ttarget_site_duplication\t24900\t24904\t.\t?\t.\tParent=repeat_region1
chr1\tLTRharvest\tLTR_retrotransposon\t25000\t27882\t.\t?\t.\tID=LTR_retrotransposon1;Parent=repeat_region1;ltr_similarity=90.87;seq_number=0
chr1\tLTRharvest\tlong_terminal_repeat\t25000\t25251\t.\t?\t.\tParent=LTR_retrotransposon1
chr1\tLTRharvest\tlong_terminal_repeat\t27638\t27882\t.\t?\t.\tParent=LTR_retrotransposon1
chr1\tLTRharvest\ttarget_site_duplication\t27883\t27887\t.\t?\t.\tParent=repeat_region1
chr2\tLTRharvest\trepeat_region\t34800\t36997\t.\t?\t.\tID=repeat_region2
chr2\tLTRharvest\tLTR_retrotransposon\t34838\t36897\t.\t?\t.\tID=LTR_retrotransposon2;Parent=repeat_region2;ltr_similarity=87.03;seq_number=1
chr2\tLTRharvest\tlong_terminal_repeat\t34838\t35022\t.\t?\t.\tParent=LTR_retrotransposon2
chr2\tLTRharvest\tlong_terminal_repeat\t36714\t36897\t.\t?\t.\tParent=LTR_retrotransposon2
"""


@pytest.fixture
def two_element_gff3(tmp_path: Path) -> Path:
    gff3 = tmp_path / "toy.gff3"
    gff3.write_text(TWO_ELEMENT_GFF3)
    return gff3


# ---------------------------------------------------------------------
# parse_ltrharvest_gff3
# ---------------------------------------------------------------------
def test_parse_reads_span_arms_similarity_and_seq_number(
    two_element_gff3: Path,
) -> None:
    """Every SCN field comes off the GFF3 in one pass, arms in left/right order."""
    records = parse_ltrharvest_gff3(two_element_gff3)
    assert len(records) == 2
    first = records[0]
    assert first.seqname == "chr1"
    assert first.seq_number == 0
    assert (first.ret_start, first.ret_end) == (25000, 27882)
    assert (first.left_start, first.left_end) == (25000, 25251)
    assert (first.right_start, first.right_end) == (27638, 27882)
    assert first.similarity == "90.87"


def test_parse_assigns_arms_to_the_enclosing_element(two_element_gff3: Path) -> None:
    """A second element's arms must not leak onto the first."""
    records = parse_ltrharvest_gff3(two_element_gff3)
    second = records[1]
    assert second.seqname == "chr2"
    assert second.seq_number == 1
    assert (second.left_start, second.right_end) == (34838, 36897)


def test_parse_ignores_repeat_region_and_tsd_features(two_element_gff3: Path) -> None:
    """Only LTR_retrotransposon and its long_terminal_repeat children matter.

    The fixture interleaves target_site_duplication rows between the element and
    its arms precisely so a naive "next two features" reader would fail here.
    """
    records = parse_ltrharvest_gff3(two_element_gff3)
    assert [r.ret_start for r in records] == [25000, 34838]


def test_parse_rejects_an_element_with_a_missing_arm(tmp_path: Path) -> None:
    """A one-armed element means a corrupt GFF3; fail loudly rather than drop it.

    Silently skipping would understate the candidate pool, and an SCN short of
    rows is indistinguishable from a genome with fewer LTR elements.
    """
    gff3 = tmp_path / "truncated.gff3"
    gff3.write_text(
        "chr1\tLTRharvest\tLTR_retrotransposon\t100\t500\t.\t?\t.\t"
        "ID=LTR_retrotransposon1;ltr_similarity=90.00;seq_number=0\n"
        "chr1\tLTRharvest\tlong_terminal_repeat\t100\t200\t.\t?\t.\t"
        "Parent=LTR_retrotransposon1\n"
    )
    with pytest.raises(ValueError, match="LTR_retrotransposon1"):
        parse_ltrharvest_gff3(gff3)


def test_parse_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="LTRharvest GFF3 not found"):
        parse_ltrharvest_gff3(tmp_path / "absent.gff3")


# ---------------------------------------------------------------------
# format_scn_row
# ---------------------------------------------------------------------
def test_format_scn_row_emits_eleven_fields_separated_by_two_spaces() -> None:
    """LTRharvest writes ``field  field`` (two spaces); lengths are computed.

    Pinned against a real row from Desmodus rotundus:
        19251  28845  9595  19251  19762  512  28330  28845  516  94.96  0
    """
    record = LtrharvestRecord(
        seqname="chr1",
        seq_number=0,
        ret_start=19251,
        ret_end=28845,
        left_start=19251,
        left_end=19762,
        right_start=28330,
        right_end=28845,
        similarity="94.96",
    )
    assert format_scn_row(record) == (
        "19251  28845  9595  19251  19762  512  28330  28845  516  94.96  0"
    )


def test_format_scn_row_preserves_similarity_text_verbatim() -> None:
    """``90.90`` must not be renormalised to ``90.9``.

    Similarity is carried as the GFF3's own string rather than parsed to float,
    so trailing zeros survive the round trip and the row stays byte-identical.
    """
    record = LtrharvestRecord(
        seqname="chr1",
        seq_number=3,
        ret_start=1,
        ret_end=10,
        left_start=1,
        left_end=3,
        right_start=8,
        right_end=10,
        similarity="90.90",
    )
    assert format_scn_row(record).split("  ")[-2] == "90.90"


def test_format_scn_row_lengths_are_closed_interval_widths() -> None:
    """l = end - start + 1 for all three spans (1-based closed intervals)."""
    record = LtrharvestRecord(
        seqname="chr1",
        seq_number=0,
        ret_start=100,
        ret_end=200,
        left_start=100,
        left_end=110,
        right_start=190,
        right_end=200,
        similarity="99.00",
    )
    fields = format_scn_row(record).split("  ")
    assert fields[2] == "101"  # retrotransposon
    assert fields[5] == "11"  # left LTR
    assert fields[8] == "11"  # right LTR


# ---------------------------------------------------------------------
# seq_number_map
# ---------------------------------------------------------------------
def test_seq_number_map_translates_scn_seq_nr_to_chromosome(
    two_element_gff3: Path,
) -> None:
    """This map replaces the ``.des`` file the prefilter used to read.

    Antrozous's ``.des`` is zero bytes, so deriving the mapping from the GFF3
    keeps the prefilter working without the suffix-array index.
    """
    records = parse_ltrharvest_gff3(two_element_gff3)
    assert seq_number_map(records) == {0: "chr1", 1: "chr2"}


def test_seq_number_map_is_empty_for_no_records() -> None:
    assert seq_number_map([]) == {}


# ---------------------------------------------------------------------
# write_scn
# ---------------------------------------------------------------------
def test_write_scn_emits_comment_header_then_rows(
    two_element_gff3: Path, tmp_path: Path
) -> None:
    """LTR_retriever skips ``#`` lines, so the header is documentation only."""
    out = tmp_path / "toy.scn"
    written = write_scn(parse_ltrharvest_gff3(two_element_gff3), out, two_element_gff3)
    assert written == 2
    lines = out.read_text().splitlines()
    assert lines[0].startswith("#")
    data = [line for line in lines if not line.startswith("#")]
    assert len(data) == 2
    assert data[0].startswith("25000  27882  2883  ")


def test_write_scn_header_documents_the_column_layout(
    two_element_gff3: Path, tmp_path: Path
) -> None:
    """The s(ret)/e(ret)/... legend is what makes a bare SCN readable."""
    out = tmp_path / "toy.scn"
    write_scn(parse_ltrharvest_gff3(two_element_gff3), out, two_element_gff3)
    assert "# s(ret) e(ret) l(ret)" in out.read_text()


def test_write_scn_creates_parent_directories(
    two_element_gff3: Path, tmp_path: Path
) -> None:
    out = tmp_path / "nested" / "deeper" / "toy.scn"
    write_scn(parse_ltrharvest_gff3(two_element_gff3), out, two_element_gff3)
    assert out.is_file()


def test_main_cli_writes_the_scn(two_element_gff3: Path, tmp_path: Path) -> None:
    out = tmp_path / "cli.scn"
    rc = scn_main(
        ["--ltrharvest-gff3", str(two_element_gff3), "--output-scn", str(out)]
    )
    assert rc == 0
    data = [line for line in out.read_text().splitlines() if not line.startswith("#")]
    assert len(data) == 2


# ---------------------------------------------------------------------
# Real-data regression (opt-in)
# ---------------------------------------------------------------------
_DEV_CACHE = os.environ.get("RETROSEEK_DEV_CACHE", "")


@pytest.mark.skipif(not _DEV_CACHE, reason="RETROSEEK_DEV_CACHE not set")
def test_rebuilt_scn_data_rows_are_byte_identical_to_real_ltrharvest_output() -> None:
    """The regression that justifies skipping the suffix-array rebuild.

    Desmodus rotundus is the one model genome with both a real LTRharvest SCN
    and its GFF3 on disk, so it is the only available ground truth. Header
    comments are regenerated (the original records absolute paths from the run
    that produced it) - the data rows are what LTR_retriever parses, and those
    must match exactly.
    """
    cache = Path(_DEV_CACHE)
    gff3 = cache / "results" / "tracks" / "ltrharvest" / "Desmodus_rotundus.gff3"
    real_scn = cache / "data" / "ltr_scn" / "Desmodus_rotundus.scn"
    if not gff3.is_file() or not real_scn.is_file():
        pytest.skip("Desmodus LTRharvest GFF3 / SCN not present in the dev cache")

    rebuilt = [format_scn_row(r) for r in parse_ltrharvest_gff3(gff3)]
    expected = [
        line
        for line in real_scn.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]
    assert rebuilt == expected
