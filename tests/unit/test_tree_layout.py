"""Unit tests for tree_layout.py (ADR-011).

Covers the taxonomy-derived cladogram, the user-Newick path, the coordinate
layout, and determinism - the last matters because the coordinates are a build
artifact that must not churn between runs.
"""

from __future__ import annotations

import pytest
import tree_layout as tl


@pytest.fixture
def taxonomy_tsv(tmp_path):
    """A miniature retroviral hierarchy: family -> subfamily -> genus."""
    p = tmp_path / "taxonomy.tsv"
    p.write_text(
        "name\tparent\trank\n"
        "Retroviridae\t\tfamily\n"
        "Orthoretrovirinae\tRetroviridae\tsubfamily\n"
        "Alpharetrovirus\tOrthoretrovirinae\tgenus\n"
        "Betaretrovirus\tOrthoretrovirinae\tgenus\n"
        "Gammaretrovirus\tOrthoretrovirinae\tgenus\n"
        "Spumaretrovirinae\tRetroviridae\tsubfamily\n"
        "Simiispumavirus\tSpumaretrovirinae\tgenus\n",
        encoding="utf-8",
    )
    return p


class TestFromTaxonomy:
    def test_every_observed_taxon_gets_exactly_one_tip(self, taxonomy_tsv):
        taxa = ["Alpharetrovirus", "Betaretrovirus", "Simiispumavirus"]
        tree = tl.from_taxonomy(taxonomy_tsv, taxa)
        tips = sorted(t.name for t in tree.get_terminals())
        assert tips == sorted(taxa)

    def test_unobserved_taxa_are_pruned_away(self, taxonomy_tsv):
        tree = tl.from_taxonomy(taxonomy_tsv, ["Alpharetrovirus", "Betaretrovirus"])
        tips = {t.name for t in tree.get_terminals()}
        assert "Simiispumavirus" not in tips  # its whole subfamily is gone

    def test_taxon_that_is_also_an_ancestor_keeps_its_own_tip(self, taxonomy_tsv):
        """Mixed-rank axis (ADR-008): loci called at `Retroviridae` need a row
        even though it is the ancestor of the genus-level calls."""
        taxa = ["Retroviridae", "Alpharetrovirus"]
        tree = tl.from_taxonomy(taxonomy_tsv, taxa)
        tips = sorted(t.name for t in tree.get_terminals())
        assert tips == ["Alpharetrovirus", "Retroviridae"]

    def test_unknown_taxa_are_dropped_not_fatal(self, taxonomy_tsv):
        tree = tl.from_taxonomy(taxonomy_tsv, ["Alpharetrovirus", "NotATaxon"])
        tips = {t.name for t in tree.get_terminals()}
        assert tips == {"Alpharetrovirus"}

    def test_no_known_taxa_returns_none(self, taxonomy_tsv):
        assert tl.from_taxonomy(taxonomy_tsv, ["NotATaxon"]) is None


class TestFromNewick:
    @pytest.fixture
    def newick(self, tmp_path):
        p = tmp_path / "species.nwk"
        p.write_text(
            "((Mus_musculus:1,Homo_sapiens:1):1,Desmodus_rotundus:2);", encoding="utf-8"
        )
        return p

    def test_tips_match_display_names_across_separator_and_case(self, newick):
        tree = tl.from_newick(newick, ["Mus musculus", "homo sapiens"])
        assert sorted(t.name for t in tree.get_terminals()) == [
            "Mus musculus",
            "homo sapiens",
        ]

    def test_species_absent_from_the_tree_are_reported_not_invented(
        self, newick, caplog
    ):
        tl.from_newick(newick, ["Mus musculus", "Myotis myotis"])
        assert "Myotis myotis" in caplog.text

    def test_zero_overlap_fails_loudly(self, newick):
        with pytest.raises(SystemExit, match="shares no tip"):
            tl.from_newick(newick, ["Gallus gallus"])


class TestLayout:
    def test_tips_get_distinct_consecutive_rows(self, taxonomy_tsv):
        tree = tl.from_taxonomy(
            taxonomy_tsv, ["Alpharetrovirus", "Betaretrovirus", "Simiispumavirus"]
        )
        _segs, tips = tl.layout(tree, align_tips=True)
        ys = sorted(t[2] for t in tips)
        assert ys == [1.0, 2.0, 3.0]

    def test_align_tips_squares_the_leaves_off(self, taxonomy_tsv):
        tree = tl.from_taxonomy(
            taxonomy_tsv, ["Alpharetrovirus", "Betaretrovirus", "Simiispumavirus"]
        )
        _segs, tips = tl.layout(tree, align_tips=True)
        assert len({t[1] for t in tips}) == 1  # all tips share one x

    def test_layout_is_deterministic(self, taxonomy_tsv):
        taxa = ["Alpharetrovirus", "Betaretrovirus", "Simiispumavirus", "Retroviridae"]
        a = tl.layout(tl.from_taxonomy(taxonomy_tsv, taxa), align_tips=True)
        b = tl.layout(tl.from_taxonomy(taxonomy_tsv, taxa), align_tips=True)
        assert a == b

    def test_segments_are_drawn_for_every_internal_node(self, taxonomy_tsv):
        tree = tl.from_taxonomy(
            taxonomy_tsv, ["Alpharetrovirus", "Betaretrovirus", "Simiispumavirus"]
        )
        segs, _tips = tl.layout(tree, align_tips=True)
        assert segs  # non-empty
        assert all(len(s) == 4 for s in segs)


class TestWrite:
    def test_headers_are_written_even_with_no_rows(self, tmp_path):
        """An unconfigured species tree still needs its files, so the DAG holds
        and the R side can render a placeholder instead of failing."""
        tl.write(tmp_path, "species", [], [])
        segs = (tmp_path / "species.tree_segments.csv").read_text(encoding="utf-8")
        tips = (tmp_path / "species.tree_tips.csv").read_text(encoding="utf-8")
        assert segs.strip() == "x,y,xend,yend"
        assert tips.strip() == "tip,x,y"


class TestSpeciesNameCanonicalization:
    """The plot generators relabel genome stems to the config `species:` display
    names before plotting, so the tree tips must carry the SAME names or the
    panel silently degrades to 'no loci matching the tree tips'. Regression
    guard for that mismatch.
    """

    def test_display_names_match_stems_across_separators(self, tmp_path):
        nwk = tmp_path / "hosts.nwk"
        nwk.write_text("(Mus_musculus:1,Homo_sapiens:1);", encoding="utf-8")
        # what the plots use, after relabel_species(stem -> display)
        display = ["Mus musculus", "Homo sapiens"]
        tree = tl.from_newick(nwk, display)
        assert sorted(t.name for t in tree.get_terminals()) == sorted(display)
