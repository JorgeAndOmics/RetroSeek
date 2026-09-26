"""Tests for solo-LTR bait construction.

The bait is what makes the whole method work without a classifier: an LTR arm
belonging to an element that hosts a catalogued ERV locus is retroviral by
construction. Three things therefore have to be exactly right, and each has a test
below: only ERV-bearing elements contribute arms, short arms are dropped, and the
GFF3 to BED coordinate shift is applied once and in the right direction.
"""

from __future__ import annotations

from pathlib import Path

import bait_builder
import pytest

from log import PipelineError

# Two elements, four arms. Element 1 hosts an ERV locus, element 2 does not.
# Element 1's right arm is deliberately short, to be dropped by the length filter.
_GFF3 = """\
##gff-version 3
chr1\tLTRharvest\tlong_terminal_repeat\t1001\t1400\t.\t+\t.\tParent=LTR_retrotransposon1;arm=L
chr1\tLTRharvest\tlong_terminal_repeat\t2001\t2100\t.\t+\t.\tParent=LTR_retrotransposon1;arm=R
chr1\tLTRharvest\tLTR_retrotransposon\t1001\t2100\t.\t+\t.\tID=LTR_retrotransposon1
chr2\tLTRharvest\tlong_terminal_repeat\t5001\t5500\t.\t-\t.\tParent=LTR_retrotransposon2;arm=L
chr2\tLTRharvest\tlong_terminal_repeat\t9001\t9500\t.\t-\t.\tParent=LTR_retrotransposon2;arm=R
"""

_LOCI = """\
id,seqname,start,end,strand,parent,source,taxon_call
L0,chr1,1200,1900,+,LTR_retrotransposon1,ltr-flanked,Betaretrovirus
"""


@pytest.fixture
def gff3(tmp_path: Path) -> Path:
    p = tmp_path / "flanking.gff3"
    p.write_text(_GFF3)
    return p


@pytest.fixture
def loci(tmp_path: Path) -> Path:
    p = tmp_path / "loci.csv"
    p.write_text(_LOCI)
    return p


# ---- which elements contribute bait ----


def test_only_erv_bearing_elements_are_read_from_the_loci_table(loci: Path) -> None:
    assert bait_builder.erv_bearing_parents(loci) == {"LTR_retrotransposon1"}


def test_a_missing_parent_column_is_an_error_not_an_empty_set(tmp_path: Path) -> None:
    """Silently returning no bait would look like "this genome has no ERVs"."""
    p = tmp_path / "no_parent.csv"
    p.write_text("id,seqname,start,end\nL0,chr1,1,2\n")
    with pytest.raises(PipelineError, match="parent"):
        bait_builder.erv_bearing_parents(p)


def test_arms_of_elements_without_an_erv_locus_are_excluded(
    gff3: Path, loci: Path
) -> None:
    parents = bait_builder.erv_bearing_parents(loci)
    arms = bait_builder.select_bait(
        bait_builder.parse_arms(gff3), parents, min_length=100
    )
    assert {a.parent for a in arms} == {"LTR_retrotransposon1"}


# ---- the length filter ----


def test_arms_shorter_than_the_minimum_are_dropped(gff3: Path, loci: Path) -> None:
    parents = bait_builder.erv_bearing_parents(loci)
    arms = bait_builder.select_bait(
        bait_builder.parse_arms(gff3), parents, min_length=300
    )
    # Element 1's L arm is 400 bp and survives; its R arm is 100 bp and does not.
    assert [a.arm for a in arms] == ["L"]


def test_the_minimum_is_inclusive(gff3: Path, loci: Path) -> None:
    parents = bait_builder.erv_bearing_parents(loci)
    arms = bait_builder.select_bait(
        bait_builder.parse_arms(gff3), parents, min_length=400
    )
    assert len(arms) == 1, "a 400 bp arm must pass a 400 bp minimum"


# ---- coordinates ----


def test_gff3_one_based_inclusive_becomes_bed_zero_based_half_open(
    gff3: Path, loci: Path, tmp_path: Path
) -> None:
    parents = bait_builder.erv_bearing_parents(loci)
    arms = bait_builder.select_bait(
        bait_builder.parse_arms(gff3), parents, min_length=300
    )
    out = tmp_path / "bait.bed"
    bait_builder.write_bed(arms, out)
    fields = out.read_text().rstrip("\n").split("\t")
    # GFF3 1001..1400 is BED 1000..1400, and the arm is 400 bp either way.
    assert fields[0] == "chr1"
    assert fields[1] == "1000"
    assert fields[2] == "1400"
    assert int(fields[2]) - int(fields[1]) == 400


def test_the_bed_name_carries_sequence_element_and_arm(
    gff3: Path, loci: Path, tmp_path: Path
) -> None:
    """The name is the only route back from a blastn hit to its seeding element."""
    parents = bait_builder.erv_bearing_parents(loci)
    arms = bait_builder.select_bait(
        bait_builder.parse_arms(gff3), parents, min_length=300
    )
    out = tmp_path / "bait.bed"
    bait_builder.write_bed(arms, out)
    assert out.read_text().split("\t")[3] == "chr1|LTR_retrotransposon1|L"


def test_bed_rows_are_six_columns_on_the_plus_strand(
    gff3: Path, loci: Path, tmp_path: Path
) -> None:
    """Bait is extracted as written (+ strand), not reverse-complemented.

    blastn searches both strands anyway, and forcing + keeps the arm's
    coordinates readable.
    """
    parents = bait_builder.erv_bearing_parents(loci)
    arms = bait_builder.select_bait(
        bait_builder.parse_arms(gff3), parents, min_length=300
    )
    out = tmp_path / "bait.bed"
    bait_builder.write_bed(arms, out)
    fields = out.read_text().rstrip("\n").split("\t")
    assert len(fields) == 6
    assert fields[5] == "+"


# ---- parsing robustness ----


def test_non_arm_features_are_ignored(gff3: Path) -> None:
    """The track carries the parent LTR_retrotransposon lines too."""
    arms = list(bait_builder.parse_arms(gff3))
    assert len(arms) == 4
    assert all(a.arm in {"L", "R"} for a in arms)


def test_an_arm_without_an_arm_attribute_is_kept_and_marked(tmp_path: Path) -> None:
    p = tmp_path / "noarm.gff3"
    p.write_text(
        "chr1\tLTRharvest\tlong_terminal_repeat\t1\t400\t.\t+\t.\tParent=LTR_retrotransposon1\n"
    )
    arms = list(bait_builder.parse_arms(p))
    assert [a.arm for a in arms] == ["?"]


def test_an_empty_gff3_yields_no_arms(tmp_path: Path) -> None:
    p = tmp_path / "empty.gff3"
    p.write_text("##gff-version 3\n")
    assert list(bait_builder.parse_arms(p)) == []
