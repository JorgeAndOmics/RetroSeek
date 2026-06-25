"""Unit tests for taxonomy_lca: the LCA primitives and GFF3 evidence parsing.

These pin the core classification contract: a single genus stays a genus, a
cross-genus set backs off to the correct higher rank, and the all-genera cloud
collapses to family level (honest abstention) rather than a forced call.
"""

from __future__ import annotations

import taxonomy_lca as tlca


class TestLca:
    def test_single_genus_returns_itself(self) -> None:
        assert tlca.lca({"Gammaretrovirus"}) == "Gammaretrovirus"

    def test_two_orthoretro_genera_back_off_to_subfamily(self) -> None:
        assert tlca.lca({"Betaretrovirus", "Gammaretrovirus"}) == "Orthoretrovirinae"

    def test_orthoretro_plus_spuma_backs_off_to_family(self) -> None:
        assert tlca.lca({"Gammaretrovirus", "Spumaretrovirus"}) == "Retroviridae"

    def test_all_seven_genera_collapse_to_family(self) -> None:
        all_genera = {
            "Alpharetrovirus",
            "Betaretrovirus",
            "Gammaretrovirus",
            "Deltaretrovirus",
            "Epsilonretrovirus",
            "Lentivirus",
            "Spumaretrovirus",
        }
        assert tlca.lca(all_genera) == "Retroviridae"

    def test_unknown_labels_ignored(self) -> None:
        assert tlca.lca({"Gammaretrovirus", "Metaviridae"}) == "Gammaretrovirus"

    def test_no_known_labels_is_unclassified(self) -> None:
        assert tlca.lca({"Metaviridae", "Bel-Pao"}) == tlca.UNCLASSIFIED


class TestWeightedLca:
    def test_strong_specific_hit_resolves_to_genus(self) -> None:
        # Gamma at 600, others far weaker -> Gammaretrovirus (not backed off).
        node, conf = tlca.weighted_lca(
            [("Gammaretrovirus", 600.0), ("Betaretrovirus", 80.0), ("Lentivirus", 75.0)]
        )
        assert node == "Gammaretrovirus"
        assert conf == 1.0

    def test_near_tie_backs_off_to_subfamily(self) -> None:
        # Gamma 600 and Beta 590 are within 10% -> both kept -> Orthoretrovirinae.
        node, conf = tlca.weighted_lca(
            [("Gammaretrovirus", 600.0), ("Betaretrovirus", 590.0)]
        )
        assert node == "Orthoretrovirinae"
        assert 0.4 < conf < 0.6

    def test_empty_is_unclassified(self) -> None:
        assert tlca.weighted_lca([]) == (tlca.UNCLASSIFIED, 0.0)

    def test_unknown_genera_ignored(self) -> None:
        node, _ = tlca.weighted_lca(
            [("Metaviridae", 900.0), ("Gammaretrovirus", 300.0)]
        )
        assert node == "Gammaretrovirus"


class TestLoadTaxonomy:
    def test_data_derived_hierarchy_drives_lca(self, tmp_path) -> None:
        # A taxonomy.tsv with a different shape than the built-in default proves
        # the hierarchy is loaded from data, not hard-coded.
        tsv = tmp_path / "taxonomy.tsv"
        tsv.write_text(
            "name\tparent\trank\n"
            "Retroviridae\t\tfamily\n"
            "Orthoretrovirinae\tRetroviridae\tsubfamily\n"
            "Gammaretrovirus\tOrthoretrovirinae\tgenus\n"
            "Betaretrovirus\tOrthoretrovirinae\tgenus\n"
        )
        saved = dict(tlca.RETRO_PARENT)
        try:
            tlca.load_taxonomy(tsv)
            assert (
                tlca.lca({"Gammaretrovirus", "Betaretrovirus"}) == "Orthoretrovirinae"
            )
            assert tlca.rank_of("Orthoretrovirinae") == "subfamily"
            # a genus absent from the loaded file is now unknown -> ignored
            assert tlca.lca({"Gammaretrovirus", "Spumaretrovirus"}) == "Gammaretrovirus"
        finally:
            tlca.RETRO_PARENT = saved  # restore built-in default for other tests


class TestRankOf:
    def test_ranks(self) -> None:
        assert tlca.rank_of("Gammaretrovirus") == "genus"
        assert tlca.rank_of("Orthoretrovirinae") == "subfamily"
        assert tlca.rank_of("Retroviridae") == "family"


class TestParseGenusSet:
    def test_decodes_escaped_semicolons(self) -> None:
        field = "Betaretrovirus%3b Deltaretrovirus%3b Gammaretrovirus"
        assert tlca.parse_genus_set(field) == {
            "Betaretrovirus",
            "Deltaretrovirus",
            "Gammaretrovirus",
        }

    def test_single_value(self) -> None:
        assert tlca.parse_genus_set("Gammaretrovirus") == {"Gammaretrovirus"}
