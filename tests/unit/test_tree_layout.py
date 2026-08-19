"""Unit tests for tree_layout.py (ADR-011).

Covers the taxonomy-derived cladogram, the user-Newick path, the coordinate
layout, and determinism - the last matters because the coordinates are a build
artifact that must not churn between runs.
"""

from __future__ import annotations

from io import StringIO
from pathlib import Path

import pytest
import tree_layout as tl
from Bio import Phylo


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


# ---------------------------------------------------------------------
# uninformative_branch_lengths (ADR-014)
#
# A cladogram is routinely supplied where a timetree is meant - the shipped
# hosts.nwk is all 1s and 2s. Topology is still usable, but a reader seeing
# branch lengths assumes they carry divergence information, so the condition is
# detected and warned about rather than passed through silently.
# ---------------------------------------------------------------------
def _tree(newick: str):
    return Phylo.read(StringIO(newick), "newick")


def test_whole_number_branch_lengths_are_uninformative() -> None:
    """The shipped hosts.nwk shape: real topology, placeholder lengths."""
    assert tl.uninformative_branch_lengths(_tree("(((A:1,M:1):1,D:2):1,(H:1,Mu:1):2);"))


def test_absent_branch_lengths_are_uninformative() -> None:
    assert tl.uninformative_branch_lengths(_tree("((A,B),(C,D));"))


def test_real_divergence_estimates_are_informative() -> None:
    """A genuine timetree carries fractional lengths and must not be flagged."""
    assert not tl.uninformative_branch_lengths(
        _tree("((A:12.4,B:12.4):43.1,(C:55.2,D:55.2):0.3);")
    )


# ---------------------------------------------------------------------
# Alias-aware tip matching
#
# Real host trees are exported from assembly pipelines and carry the genome
# DIRECTORY name, which in a real study bears no textual relation to the
# display name: `GCF_014176215.1_mMyoMyo1` is "Myotis myotis". Folding
# underscores and case cannot bridge that, but the config `species:` block
# already states the correspondence, so matching consults it.
# ---------------------------------------------------------------------
SPECIES_MAP = {
    "GCF_014176215.1_mMyoMyo1": "Myotis myotis",
    "Antrozous_pallidus": "Antrozous pallidus",
    "GCA_004026805.1_ASM402680v1": "Desmodus rotundus",
}


def test_alias_index_maps_stems_and_display_names_to_the_display_name() -> None:
    """Both spellings resolve to the display name the catalog uses."""
    idx = tl.build_alias_index(SPECIES_MAP)
    assert idx[tl._normalize("GCF_014176215.1_mMyoMyo1")] == "Myotis myotis"
    assert idx[tl._normalize("Myotis myotis")] == "Myotis myotis"
    assert idx[tl._normalize("myotis_myotis")] == "Myotis myotis"


def test_alias_index_can_canonicalise_to_the_stem_instead() -> None:
    """The co-phylogeny compares against trees whose tips ARE genome stems,
    so it needs the mapping pointing the other way."""
    idx = tl.build_alias_index(SPECIES_MAP, canonical="stem")
    assert idx[tl._normalize("Myotis myotis")] == "GCF_014176215.1_mMyoMyo1"
    assert idx[tl._normalize("GCF_014176215.1_mMyoMyo1")] == "GCF_014176215.1_mMyoMyo1"


def test_from_newick_matches_a_tree_labelled_with_genome_directories(
    tmp_path: Path,
) -> None:
    """The real-workflow case: an assembly-derived tree, arbitrary stem labels.

    Without the alias map every tip fails to match and the panel silently
    renders "no loci matching the tree tips" - a blank figure, not an error.
    """
    nwk = tmp_path / "assembly.tre"
    nwk.write_text(
        "((GCF_014176215.1_mMyoMyo1:12.4,Antrozous_pallidus:12.4):8.1,"
        "GCA_004026805.1_ASM402680v1:20.5);"
    )
    tree = tl.from_newick(
        nwk,
        ["Myotis myotis", "Antrozous pallidus", "Desmodus rotundus"],
        aliases=tl.build_alias_index(SPECIES_MAP),
    )
    # tips are renamed to the display names the catalog keys on
    assert sorted(x.name for x in tree.get_terminals()) == [
        "Antrozous pallidus",
        "Desmodus rotundus",
        "Myotis myotis",
    ]


def test_from_newick_still_matches_display_names_without_aliases(
    tmp_path: Path,
) -> None:
    """The existing behaviour must not regress when no map is supplied."""
    nwk = tmp_path / "plain.nwk"
    nwk.write_text("(Antrozous_pallidus:1,Mus_musculus:1);")
    tree = tl.from_newick(nwk, ["Antrozous pallidus", "Mus musculus"])
    assert sorted(x.name for x in tree.get_terminals()) == [
        "Antrozous pallidus",
        "Mus musculus",
    ]


def test_from_newick_prunes_tips_the_study_does_not_include(tmp_path: Path) -> None:
    """A 100-species reference timetree pruned down to the study's genomes."""
    nwk = tmp_path / "big.tre"
    nwk.write_text(
        "((GCF_014176215.1_mMyoMyo1:1,Bos_taurus:1):1,(Gallus_gallus:1,"
        "Antrozous_pallidus:1):1);"
    )
    tree = tl.from_newick(
        nwk,
        ["Myotis myotis", "Antrozous pallidus"],
        aliases=tl.build_alias_index(SPECIES_MAP),
    )
    assert sorted(x.name for x in tree.get_terminals()) == [
        "Antrozous pallidus",
        "Myotis myotis",
    ]
