"""Tests for cutting the solo-LTR tree into families.

The tree below is small enough to work out by hand:

    root
    |-- 0.5 -- A: (FLANK a:0.05, FLANK a:0.05)         diameter 0.10
    `-- 0.3 -- BC
               |-- 0.4 -- B: (SOLO 1:0.02, SOLO 2:0.03)   diameter 0.05
               `-- 0.1 -- C: (SOLO 3:0.05, FLANK c:0.05)  diameter 0.10

BC's diameter is the path SOLO 2 -> BC -> SOLO 3: 0.03 + 0.4 + 0.1 + 0.05 = 0.58,
so with a 0.2 cut the families are exactly A, B and C.
"""

from __future__ import annotations

from io import StringIO

import pytest
import tree_families

pytest.importorskip("Bio")
from Bio import Phylo

_NEWICK = (
    "((FLANK__a_e1_L:0.05,FLANK__a_e1_R:0.05):0.5,"
    "((SOLO__b_1:0.02,SOLO__b_2:0.03):0.4,(SOLO__c_3:0.05,FLANK__c_e3_L:0.05):0.1):0.3);"
)


def _tree():
    return Phylo.read(StringIO(_NEWICK), "newick")


def test_diameter_is_the_longest_path_inside_a_clade() -> None:
    tree = _tree()
    tree_families.annotate_diameters(tree)
    bc = tree.common_ancestor("SOLO__b_2", "SOLO__c_3")
    assert bc.diameter == pytest.approx(0.58)


def test_the_cut_returns_maximal_clades_within_the_distance() -> None:
    families = tree_families.cut_families(_tree(), max_distance=0.2)
    members = [sorted(t.name for t in f.clade.get_terminals()) for f in families]
    assert sorted(members) == [
        ["FLANK__a_e1_L", "FLANK__a_e1_R"],
        ["FLANK__c_e3_L", "SOLO__c_3"],
        ["SOLO__b_1", "SOLO__b_2"],
    ]


def test_a_large_enough_distance_makes_the_whole_tree_one_family() -> None:
    families = tree_families.cut_families(_tree(), max_distance=10.0)
    assert len(families) == 1
    assert families[0].n_tips == 6


def test_each_family_gets_the_right_kind() -> None:
    families = {
        tuple(sorted(t.name for t in f.clade.get_terminals())): f.kind
        for f in tree_families.cut_families(_tree(), max_distance=0.2)
    }
    assert families[("FLANK__a_e1_L", "FLANK__a_e1_R")] == tree_families.NO_SOLO
    assert families[("SOLO__b_1", "SOLO__b_2")] == tree_families.NO_INTACT
    assert families[("FLANK__c_e3_L", "SOLO__c_3")] == tree_families.WITH_INTACT


def test_family_ids_are_stable_and_start_with_the_largest() -> None:
    first = [f.family_id for f in tree_families.cut_families(_tree(), 0.2)]
    second = [f.family_id for f in tree_families.cut_families(_tree(), 0.2)]
    assert first == second == ["F001", "F002", "F003"]


def test_every_tip_is_assigned_to_exactly_one_family() -> None:
    families = tree_families.cut_families(_tree(), 0.2)
    mapping = tree_families.tip_families(families)
    assert len(mapping) == 6


def test_the_solo_only_tree_keeps_solos_and_their_distances() -> None:
    """Pruning must not change the path between the remaining solos."""
    tree = _tree()
    before = tree.distance("SOLO__b_1", "SOLO__c_3")
    pruned = tree_families.solo_only_tree(tree)
    assert pruned is not None
    assert sorted(t.name for t in pruned.get_terminals()) == [
        "SOLO__b_1",
        "SOLO__b_2",
        "SOLO__c_3",
    ]
    assert pruned.distance("SOLO__b_1", "SOLO__c_3") == pytest.approx(before)


def test_the_solo_only_tree_is_none_when_too_few_solos() -> None:
    tree = Phylo.read(
        StringIO("((FLANK__a_L:1,FLANK__a_R:1):1,SOLO__x_1:1);"), "newick"
    )
    assert tree_families.solo_only_tree(tree) is None


def test_showcase_picks_both_kinds_and_skips_tiny_families() -> None:
    """Families under three tips are not worth drawing as a tree of their own."""
    newick = (
        "(((SOLO__a_1:0.01,SOLO__a_2:0.01):0.01,SOLO__a_3:0.01):1,"
        "((SOLO__b_4:0.01,FLANK__b_e_L:0.01):0.01,FLANK__b_e_R:0.01):1,"
        "(SOLO__c_5:0.01,SOLO__c_6:0.01):1);"
    )
    families = tree_families.cut_families(Phylo.read(StringIO(newick), "newick"), 0.2)
    chosen = tree_families.showcase(families, per_kind=2)
    assert sorted(f.kind for f in chosen) == [
        tree_families.NO_INTACT,
        tree_families.WITH_INTACT,
    ]
    assert all(f.n_tips >= 3 for f in chosen)
