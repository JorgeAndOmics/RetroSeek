"""Tests for the solo-LTR tree views.

Every figure page about the tree is drawn from these files, so they must agree with
each other: each tip carries the class its name encodes (or the colours lie), every
tip belongs to one family, and the solo-only tree holds exactly the solos.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest
import solo_tree_layout

pytest.importorskip("Bio")

# Two families at a 0.2 cut: flanking arms a, and a mixed family b of three solos
# and one flanking arm. Branch lengths vary so x is not trivially uniform.
_NEWICK = (
    "((FLANK__c_e1_L:0.05,FLANK__c_e1_R:0.06):0.9,"
    "(((SOLO__c_10:0.02,SOLO__c_11:0.03):0.02,SOLO__c_12:0.04):0.01,FLANK__c_e2_L:0.03):0.8);"
)


def _read(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


@pytest.fixture
def views(tmp_path: Path) -> Path:
    tree = tmp_path / "t.treefile"
    tree.write_text(_NEWICK)
    prefix = tmp_path / "G"
    solo_tree_layout.write_views(
        tree, prefix, tmp_path / "G.solos.treefile", max_distance=0.2, per_kind=3
    )
    return tmp_path


def test_every_tip_is_written_with_its_class_and_family(views: Path) -> None:
    tips = _read(views / "G.tree_tips.csv")
    assert sorted(t["class"] for t in tips) == [
        "FLANK",
        "FLANK",
        "FLANK",
        "SOLO",
        "SOLO",
        "SOLO",
    ]
    assert all(t["family"].startswith("F") for t in tips)


def test_tips_occupy_one_row_each(views: Path) -> None:
    """y is the row a tip is drawn on; two tips on one row would overprint."""
    tips = _read(views / "G.tree_tips.csv")
    assert sorted(float(t["y"]) for t in tips) == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]


def test_branch_lengths_are_kept(views: Path) -> None:
    """On an LTR tree divergence is the signal, so tips must not be squared off."""
    tips = _read(views / "G.tree_tips.csv")
    assert len({float(t["x"]) for t in tips}) > 1


def test_the_solo_only_view_holds_exactly_the_solos(views: Path) -> None:
    tips = _read(views / "G.solo_tree_tips.csv")
    assert sorted(t["tip"] for t in tips) == ["SOLO__c_10", "SOLO__c_11", "SOLO__c_12"]
    assert (views / "G.solos.treefile").read_text().strip().endswith(";")


def test_the_family_table_names_each_kind(views: Path) -> None:
    kinds = {row["kind"] for row in _read(views / "G.tree_families.csv")}
    assert kinds == {"no_solo", "with_intact"}


def test_only_showcased_families_get_their_own_subtree(views: Path) -> None:
    families = _read(views / "G.tree_families.csv")
    shown = {row["family"] for row in families if row["shown"] == "True"}
    drawn = {row["family"] for row in _read(views / "G.family_tree_tips.csv")}
    assert drawn == shown
    # The flanking-only family has no solos, so it is never showcased.
    assert all(row["kind"] != "no_solo" for row in families if row["family"] in shown)


def test_too_few_solos_still_writes_every_file(tmp_path: Path) -> None:
    """The figure must show a placeholder, not fail, when a view is empty."""
    tree = tmp_path / "t.treefile"
    tree.write_text("((FLANK__a_L:0.1,FLANK__a_R:0.1):0.5,SOLO__x_1:0.3);")
    prefix = tmp_path / "G"
    solo_tree_layout.write_views(tree, prefix, tmp_path / "G.solos.treefile", 0.2, 3)
    assert _read(tmp_path / "G.solo_tree_tips.csv") == []
    assert (tmp_path / "G.solos.treefile").exists()
