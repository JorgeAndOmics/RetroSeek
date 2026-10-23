"""Tests for the upstream Dfam filter of LTR elements (ltr_arm_filter.py, ADR-026).

An LTRdigest element whose LTR arms are really LINE, SINE or DNA repeats (two
copies of an L1 3' end, two Alus) is not an ERV, yet an ERV probe hit between them
makes it look like one and its arms then seed the solo-LTR search. The filter drops
such elements before the loci exist. What can go wrong quietly: dropping an element
because an arm has no label (unlabelled arms must stay), keeping an element whose
other arm is a LINE, a hit touching an element by one base being missed, minus
strand hits (start after end) being ignored, and arm coordinates being off by one.
"""

from __future__ import annotations

from pathlib import Path

import ltr_arm_filter as f
import pytest
from bait_builder import Arm

from log import PipelineError

# Two elements on chr1: E1 (arms 101-200 and 801-900) and E2 (arms 3001-3100 and
# 3801-3900); one on chr2. Element rows carry the container's id as Parent, as in
# LTRdigest's output; arms and domains carry the element's id.
GFF3 = """\
##gff-version 3
chr1\tLTRharvest\trepeat_region\t96\t905\t.\t?\t.\tID=repeat_region1
chr1\tLTRharvest\ttarget_site_duplication\t96\t100\t.\t?\t.\tParent=repeat_region1
chr1\tLTRharvest\tLTR_retrotransposon\t101\t900\t.\t+\t.\tID=LTR_retrotransposon1;Parent=repeat_region1
chr1\tLTRharvest\tlong_terminal_repeat\t801\t900\t.\t+\t.\tParent=LTR_retrotransposon1
chr1\tLTRharvest\tlong_terminal_repeat\t101\t200\t.\t+\t.\tParent=LTR_retrotransposon1
chr1\tLTRdigest\tprotein_match\t300\t500\t.\t+\t.\tParent=LTR_retrotransposon1;name=RVT_1
chr1\tLTRharvest\tLTR_retrotransposon\t3001\t3900\t.\t-\t.\tID=LTR_retrotransposon2;Parent=repeat_region2
chr1\tLTRharvest\tlong_terminal_repeat\t3001\t3100\t.\t-\t.\tParent=LTR_retrotransposon2
chr1\tLTRharvest\tlong_terminal_repeat\t3801\t3900\t.\t-\t.\tParent=LTR_retrotransposon2
chr2\tLTRharvest\tLTR_retrotransposon\t11\t60\t.\t+\t.\tID=LTR_retrotransposon3;Parent=repeat_region3
chr2\tLTRharvest\tlong_terminal_repeat\t11\t20\t.\t+\t.\tParent=LTR_retrotransposon3
chr2\tLTRharvest\tlong_terminal_repeat\t51\t60\t.\t+\t.\tParent=LTR_retrotransposon3
"""


@pytest.fixture
def gff3(tmp_path: Path) -> Path:
    path = tmp_path / "ltrdigest.gff3"
    path.write_text(GFF3)
    return path


# ---------------------------------------------------------------------------
# The decision: which classes are not LTRs, and when an element goes
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "label", ["LTR/ERVK", "LTR/ERV1", "LTR/Gypsy", "LTR", "ltr/ervl"]
)
def test_every_ltr_class_is_kept(label):
    assert not f.is_non_ltr(label)


@pytest.mark.parametrize("label", ["", "Unknown", "unknown", None])
def test_an_unlabelled_arm_is_not_a_reason_to_drop(label):
    assert not f.is_non_ltr(label)


@pytest.mark.parametrize(
    "label", ["LINE", "LINE/L1", "SINE", "SINE/Alu", "DNA", "DNA/hAT", "Satellite"]
)
def test_line_sine_dna_and_other_non_ltr_classes_are_dropped(label):
    assert f.is_non_ltr(label)


def test_an_element_with_two_ltr_arms_stays():
    assert f.decide(["LTR/ERVK", "LTR/ERVK"]) == (False, "")


def test_either_arm_being_non_ltr_drops_the_element():
    drop, reason = f.decide(["LTR/ERV1", "LINE/L1"])
    assert drop
    assert "LINE/L1" in reason


def test_a_non_ltr_arm_beside_an_unlabelled_arm_still_drops():
    assert f.decide(["", "SINE/Alu"])[0]


def test_an_unlabelled_arm_beside_an_ltr_arm_stays():
    assert f.decide(["", "LTR/ERVK"]) == (False, "")


def test_two_unlabelled_arms_stay():
    assert f.decide(["", ""]) == (False, "")


# ---------------------------------------------------------------------------
# Reading the LTRdigest GFF3
# ---------------------------------------------------------------------------
def test_each_element_gets_its_two_arms_in_genomic_order(gff3):
    arms = [a for a in f.element_arms(gff3) if a.parent == "LTR_retrotransposon1"]
    assert [(a.start, a.end, a.arm) for a in arms] == [(101, 200, "1"), (801, 900, "2")]


def test_arm_names_are_unique_and_name_the_element(gff3):
    names = [a.name for a in f.element_arms(gff3)]
    assert len(names) == len(set(names)) == 6
    assert "chr2|LTR_retrotransposon3|1" in names


def test_domains_and_containers_are_not_arms(gff3):
    assert {a.parent for a in f.element_arms(gff3)} == {
        "LTR_retrotransposon1",
        "LTR_retrotransposon2",
        "LTR_retrotransposon3",
    }


def test_element_spans_come_from_the_element_row_not_the_container(gff3):
    spans = f.element_spans(gff3)
    assert spans["LTR_retrotransposon1"] == ("chr1", 101, 900)
    assert len(spans) == 3


# ---------------------------------------------------------------------------
# Which elements bear a hit (a superset of the R side's strand-aware overlap)
# ---------------------------------------------------------------------------
def test_a_hit_inside_an_element_marks_it(gff3):
    spans = f.element_spans(gff3)
    assert f.overlapping_elements(spans, [("chr1", 300, 400)]) == {
        "LTR_retrotransposon1"
    }


def test_a_hit_touching_an_element_by_one_base_marks_it():
    spans = {"E": ("chr1", 101, 900)}
    assert f.overlapping_elements(spans, [("chr1", 50, 101)]) == {"E"}
    assert f.overlapping_elements(spans, [("chr1", 900, 950)]) == {"E"}


def test_a_hit_one_base_outside_does_not_mark_it():
    spans = {"E": ("chr1", 101, 900)}
    assert f.overlapping_elements(spans, [("chr1", 50, 100)]) == set()
    assert f.overlapping_elements(spans, [("chr1", 901, 950)]) == set()


def test_a_hit_on_another_sequence_does_not_mark_it():
    assert (
        f.overlapping_elements({"E": ("chr1", 101, 900)}, [("chr2", 300, 400)]) == set()
    )


def test_minus_strand_hits_given_end_first_still_count():
    assert f.overlapping_elements({"E": ("chr1", 101, 900)}, [("chr1", 400, 300)]) == {
        "E"
    }


def test_a_hit_spanning_several_elements_marks_all_of_them(gff3):
    spans = f.element_spans(gff3)
    assert f.overlapping_elements(spans, [("chr1", 500, 3500)]) == {
        "LTR_retrotransposon1",
        "LTR_retrotransposon2",
    }


def test_no_hits_marks_nothing(gff3):
    assert f.overlapping_elements(f.element_spans(gff3), []) == set()


# ---------------------------------------------------------------------------
# Cutting the arms out of the genome
# ---------------------------------------------------------------------------
def test_arms_are_cut_one_based_inclusive(tmp_path):
    fasta = tmp_path / "g.fa"
    fasta.write_text(
        ">chr1 some description\nACGTACGTAC\nGGTTAACCGG\n>chr2\nTTTTTAAAAA\n"
    )
    arms = [
        Arm("chr1", 1, 4, "E1", "1"),
        Arm("chr1", 9, 12, "E1", "2"),
        Arm("chr2", 6, 10, "E2", "1"),
    ]
    cut = f.extract_arms(fasta, arms)
    assert cut["chr1|E1|1"] == "ACGT"
    assert cut["chr1|E1|2"] == "ACGG"  # crosses the line break
    assert cut["chr2|E2|1"] == "AAAAA"


def test_an_arm_on_a_missing_sequence_is_an_error(tmp_path):
    fasta = tmp_path / "g.fa"
    fasta.write_text(">chr1\nACGT\n")
    with pytest.raises(PipelineError, match="chrX"):
        f.extract_arms(fasta, [Arm("chrX", 1, 3, "E", "1")])


# ---------------------------------------------------------------------------
# The verdict table
# ---------------------------------------------------------------------------
ARMS = [
    Arm("chr1", 101, 200, "E1", "1"),
    Arm("chr1", 801, 900, "E1", "2"),
    Arm("chr1", 3001, 3100, "E2", "1"),
    Arm("chr1", 3801, 3900, "E2", "2"),
]
FAMILY_OF = {a.name: "famA" for a in ARMS[:2]} | {a.name: "famB" for a in ARMS[2:]}


def test_each_element_gets_one_row_with_both_arms_labelled():
    rows = f.verdict_rows(ARMS, FAMILY_OF, {"famA": "LTR/ERVK", "famB": "LINE/L1"})
    by_element = {r["element"]: r for r in rows}
    assert len(rows) == 2
    assert by_element["E1"]["drop"] is False
    assert by_element["E2"]["drop"] is True
    assert by_element["E2"]["arm1_class"] == by_element["E2"]["arm2_class"] == "LINE/L1"
    assert by_element["E2"]["arm1_family"] == "famB"


def test_an_arm_in_no_family_is_unlabelled_and_cannot_drop_its_element():
    family_of = {ARMS[0].name: "famA"}  # the other arm was too short to be clustered
    rows = f.verdict_rows(ARMS[:2], family_of, {"famA": "LTR/ERV1"})
    assert rows[0]["arm2_family"] == ""
    assert rows[0]["arm2_class"] == ""
    assert rows[0]["drop"] is False


def test_a_family_with_no_dfam_match_is_unlabelled():
    rows = f.verdict_rows(ARMS[:2], FAMILY_OF, {})  # nothing matched any model
    assert rows[0]["drop"] is False
    assert rows[0]["arm1_class"] == ""


def test_rows_follow_the_input_element_order():
    rows = f.verdict_rows(ARMS, FAMILY_OF, {"famA": "LTR/ERVK", "famB": "LTR/ERVK"})
    assert [r["element"] for r in rows] == ["E1", "E2"]
    assert all(r["seqname"] == "chr1" for r in rows)
