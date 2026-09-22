"""Tests for the solo-LTR tree layout.

The layout is what the figure is drawn from, so two things must hold: every tip
carries the class its name encodes (or the colours lie), and the coordinates are
the same ones taxonomy/tree_layout.py produces (the code is reused, not copied).
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest
import solo_tree_layout

pytest.importorskip("Bio")

# ((F,F),(S,M)) with branch lengths, so x is not trivially uniform.
_NEWICK = (
    "((FLANK__c_e1_L:0.1,FLANK__c_e1_R:0.1):0.2,(SOLO__c_10:0.3,MONO__c_20:0.4):0.1);"
)


def _layout(tmp_path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    tree = tmp_path / "t.treefile"
    tree.write_text(_NEWICK)
    tips_csv, segs_csv = tmp_path / "tips.csv", tmp_path / "segs.csv"
    solo_tree_layout.write_layout(tree, tips_csv, segs_csv)
    with tips_csv.open() as h:
        tips = list(csv.DictReader(h))
    with segs_csv.open() as h:
        segs = list(csv.DictReader(h))
    return tips, segs


def test_every_tip_is_written_with_its_class(tmp_path: Path) -> None:
    tips, _ = _layout(tmp_path)
    assert sorted(t["class"] for t in tips) == ["FLANK", "FLANK", "MONO", "SOLO"]


def test_tips_occupy_one_row_each(tmp_path: Path) -> None:
    """y is the row a tip is drawn on; two tips on one row would overprint."""
    tips, _ = _layout(tmp_path)
    assert sorted(float(t["y"]) for t in tips) == [1.0, 2.0, 3.0, 4.0]


def test_branch_lengths_are_kept(tmp_path: Path) -> None:
    """On an LTR tree divergence is the signal, so tips must not be squared off."""
    tips, _ = _layout(tmp_path)
    assert len({float(t["x"]) for t in tips}) > 1


def test_segments_connect_every_internal_node(tmp_path: Path) -> None:
    """Three internal nodes (root plus two cherries), each one vertical spine and
    two horizontal arms."""
    _, segs = _layout(tmp_path)
    assert len(segs) == 9
