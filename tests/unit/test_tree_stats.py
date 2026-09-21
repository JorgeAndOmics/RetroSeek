"""Tests for the solo-LTR tree statistics.

These numbers are the evidence that the three LTR fates are real classes, so each
is pinned on a hand-built tree small enough to verify by eye. The Newick strings
below are drawn in comments: if a test fails, the expected answer is countable from
the drawing rather than trusted.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import tree_stats

pytest.importorskip("Bio")


def _tree(newick: str, tmp_path: Path) -> Path:
    p = tmp_path / "t.treefile"
    p.write_text(newick)
    return p


# ---- class parsing ----


def test_the_class_is_the_prefix_before_the_double_underscore() -> None:
    assert tree_stats.tip_class("FLANK__chr1_LTR_retrotransposon1_L") == "FLANK"
    assert tree_stats.tip_class("SOLO__chr1_1000") == "SOLO"
    assert tree_stats.tip_class("MONO__chr1_1000") == "MONO"


# ---- the positive control ----


def test_both_arms_of_one_element_as_sisters_are_counted(tmp_path: Path) -> None:
    """((L,R),other) - the two arms are sisters, so 1 of 1."""
    newick = "((FLANK__chr1_elemA_L,FLANK__chr1_elemA_R),SOLO__chr1_500);"
    summary = tree_stats.summarise(_tree(newick, tmp_path), permutations=2, seed=1)
    assert summary["elements_with_both_arms"] == 1
    assert summary["arms_recovered_as_sisters"] == 1
    assert summary["arm_sisterhood_fraction"] == 1.0


def test_arms_pulled_apart_are_counted_as_a_control_failure(tmp_path: Path) -> None:
    """((L,solo),R) - the arms are not sisters, so 0 of 1."""
    newick = "((FLANK__chr1_elemA_L,SOLO__chr1_500),FLANK__chr1_elemA_R);"
    summary = tree_stats.summarise(_tree(newick, tmp_path), permutations=2, seed=1)
    assert summary["elements_with_both_arms"] == 1
    assert summary["arms_recovered_as_sisters"] == 0


def test_an_element_with_only_one_arm_on_the_tree_is_not_counted(
    tmp_path: Path,
) -> None:
    """The control asks about pairs; a single arm cannot pass or fail it."""
    newick = "(FLANK__chr1_elemA_L,SOLO__chr1_500);"
    summary = tree_stats.summarise(_tree(newick, tmp_path), permutations=2, seed=1)
    assert summary["elements_with_both_arms"] == 0
    assert summary["arm_sisterhood_fraction"] == ""


# ---- clustering ----


def test_perfectly_clustered_classes_give_a_fraction_of_one(tmp_path: Path) -> None:
    """((S,S),(F,F)) - every tip's sister is its own class."""
    newick = "((SOLO__a_1,SOLO__a_2),(FLANK__b_x_L,FLANK__b_x_R));"
    summary = tree_stats.summarise(_tree(newick, tmp_path), permutations=2, seed=1)
    assert summary["same_class_sister_observed"] == 1.0


def test_perfectly_interleaved_classes_give_zero(tmp_path: Path) -> None:
    """((S,F),(S,F)) - no tip's immediate sister shares its class.

    The outer pairs are subtrees, so each tip's sister group is its partner only.
    """
    newick = "((SOLO__a_1,FLANK__b_x_L),(SOLO__a_2,FLANK__b_y_L));"
    summary = tree_stats.summarise(_tree(newick, tmp_path), permutations=2, seed=1)
    assert summary["same_class_sister_observed"] == 0.0


def test_the_null_is_reported_alongside_the_observation(tmp_path: Path) -> None:
    """A clustering number without its null is uninterpretable, so both are kept."""
    newick = "((SOLO__a_1,SOLO__a_2),(FLANK__b_x_L,FLANK__b_x_R));"
    summary = tree_stats.summarise(_tree(newick, tmp_path), permutations=5, seed=7)
    assert "same_class_sister_null_mean" in summary
    assert summary["permutations"] == 5
    assert summary["seed"] == 7


def test_the_null_is_reproducible_for_a_given_seed(tmp_path: Path) -> None:
    newick = "((SOLO__a_1,SOLO__a_2),(FLANK__b_x_L,MONO__c_9));"
    path = _tree(newick, tmp_path)
    first = tree_stats.summarise(path, permutations=10, seed=42)
    second = tree_stats.summarise(path, permutations=10, seed=42)
    assert first["same_class_sister_null_mean"] == second["same_class_sister_null_mean"]


# ---- census and adjacency ----


def test_the_tip_census_counts_each_class(tmp_path: Path) -> None:
    newick = "((SOLO__a_1,SOLO__a_2),(FLANK__b_x_L,MONO__c_9));"
    summary = tree_stats.summarise(_tree(newick, tmp_path), permutations=2, seed=1)
    assert (
        summary["n_tips"],
        summary["n_solo"],
        summary["n_flank"],
        summary["n_mono"],
    ) == (
        4,
        2,
        1,
        1,
    )


def test_adjacency_covers_every_class_pair(tmp_path: Path) -> None:
    """A fixed 3x3 shape, so the plot never has to handle a missing cell."""
    newick = "((SOLO__a_1,SOLO__a_2),(FLANK__b_x_L,MONO__c_9));"
    summary = tree_stats.summarise(_tree(newick, tmp_path), permutations=2, seed=1)
    pairs = {(row["tip_class"], row["sister_class"]) for row in summary["adjacency"]}
    assert len(pairs) == 9


def test_adjacency_enrichment_is_relative_to_class_abundance(tmp_path: Path) -> None:
    """Two solos that are each other's sister, in a tree half solo by abundance.

    Observed SOLO-beside-SOLO is 1.0, expected is 2/4, so enrichment is 2.0.
    """
    newick = "((SOLO__a_1,SOLO__a_2),(FLANK__b_x_L,MONO__c_9));"
    summary = tree_stats.summarise(_tree(newick, tmp_path), permutations=2, seed=1)
    cell = next(
        row
        for row in summary["adjacency"]
        if row["tip_class"] == "SOLO" and row["sister_class"] == "SOLO"
    )
    assert cell["observed_fraction"] == 1.0
    assert cell["expected_fraction"] == 0.5
    assert cell["enrichment"] == 2.0


# ---- output ----


def test_both_csvs_are_written_and_readable(tmp_path: Path) -> None:
    newick = "((SOLO__a_1,SOLO__a_2),(FLANK__b_x_L,FLANK__b_x_R));"
    summary = tree_stats.summarise(_tree(newick, tmp_path), permutations=3, seed=1)
    summary_csv = tmp_path / "summary.csv"
    adjacency_csv = tmp_path / "adjacency.csv"
    tree_stats.write_csvs(summary, summary_csv, adjacency_csv)

    text = summary_csv.read_text()
    assert text.startswith("metric,value")
    assert "arm_sisterhood_fraction" in text
    assert adjacency_csv.read_text().count("\n") == 10  # header plus 3x3
