"""Unit tests for the per-locus classifier internals (taxonomy_classify_loci).

These pin the data-derived, gene-agnostic contract: diagnostic genes are
discovered from the reference (never hard-coded), loci are grouped by their
LTR-element ``Parent``, and the assembled record carries the taxon-founded
structural metrics (completeness, canonical order, mosaic). A locus is 'resolved'
when its call lands on a declared axis taxon (ADR-008, rank-agnostic). A separate
test covers the gappa-output parser used by the placement branch.
"""

from __future__ import annotations

import taxonomy_classify_loci as tcl
import taxonomy_lca as tlca
import taxonomy_placement as tplace

# Default (retroviral-genus) axis used by most tests — matches the fallback taxonomy.
_AXIS = {
    "Alpharetrovirus",
    "Betaretrovirus",
    "Gammaretrovirus",
    "Deltaretrovirus",
    "Epsilonretrovirus",
    "Lentivirus",
    "Spumaretrovirus",
}


class TestAutoDiagnostic:
    def test_single_taxon_gene_is_diagnostic(self) -> None:
        gene_of = {"a1": "REX", "a2": "REX", "a3": "POL"}
        taxon_of = {
            "a1": "Deltaretrovirus",
            "a2": "Deltaretrovirus",
            "a3": "Gammaretrovirus",
        }
        diag = tcl.auto_diagnostic(gene_of, taxon_of)
        assert diag["REX"] == "Deltaretrovirus"

    def test_multi_taxon_gene_not_diagnostic(self) -> None:
        gene_of = {"a1": "POL", "a2": "POL"}
        taxon_of = {"a1": "Gammaretrovirus", "a2": "Betaretrovirus"}
        assert "POL" not in tcl.auto_diagnostic(gene_of, taxon_of)

    def test_other_catchall_never_diagnostic(self) -> None:
        gene_of = {"a1": "OTHER"}
        taxon_of = {"a1": "Lentivirus"}
        assert "OTHER" not in tcl.auto_diagnostic(gene_of, taxon_of)


class TestRefMaps:
    def test_returns_taxon_gene_maps_and_axis(self, tmp_path) -> None:
        csv = tmp_path / "retro_reference.csv"
        csv.write_text(
            "accession,taxon,gene,defline\n"
            "A1,Lentivirus,POL,pol protein\n"
            "A2,Gammaretrovirus,GAG,gag protein\n"
            "A3,Lentivirus,ENV,env protein\n",
            encoding="utf-8",
        )
        taxon_of, gene_of, axis = tcl._ref_maps(csv)
        assert taxon_of["A1"] == "Lentivirus"
        assert gene_of["A2"] == "GAG"
        # axis = the distinct declared reference taxa (any rank)
        assert axis == {"Lentivirus", "Gammaretrovirus"}


class TestBuildLoci:
    def test_features_group_by_parent(self) -> None:
        feats = [
            {
                "seqname": "chr1",
                "start": "100",
                "end": "200",
                "strand": "+",
                "gene": "POL",
                "parent": "retro1",
                "label": "MLV",
            },
            {
                "seqname": "chr1",
                "start": "210",
                "end": "300",
                "strand": "+",
                "gene": "GAG",
                "parent": "retro1",
                "label": "MLV",
            },
        ]
        loci = tcl.build_loci(feats)
        assert len(loci) == 1
        locus = loci[0]
        assert set(locus["genes"]) == {"POL", "GAG"}
        assert locus["start"] == 100
        assert locus["end"] == 300

    def test_parentless_features_become_own_locus(self) -> None:
        feats = [
            {
                "seqname": "chr1",
                "start": "1",
                "end": "9",
                "strand": "+",
                "gene": "POL",
                "parent": "",
                "label": "",
            },
            {
                "seqname": "chr1",
                "start": "20",
                "end": "29",
                "strand": "+",
                "gene": "POL",
                "parent": "",
                "label": "",
            },
        ]
        assert len(tcl.build_loci(feats)) == 2

    def test_same_gene_coords_merged_to_span(self) -> None:
        feats = [
            {
                "seqname": "chr1",
                "start": "100",
                "end": "150",
                "strand": "+",
                "gene": "POL",
                "parent": "r",
                "label": "",
            },
            {
                "seqname": "chr1",
                "start": "180",
                "end": "260",
                "strand": "+",
                "gene": "POL",
                "parent": "r",
                "label": "",
            },
        ]
        locus = tcl.build_loci(feats)[0]
        assert locus["genes"]["POL"] == (100, 260)


class TestAssembleStructure:
    def _locus(self, genes: dict[str, tuple[int, int]]) -> dict:
        return {
            "id": "L0",
            "seqname": "chr1",
            "parent": "r",
            "strand": "+",
            "start": min(s for s, _ in genes.values()),
            "end": max(e for _, e in genes.values()),
            "genes": genes,
            "probe_label_set": "",
        }

    def test_completeness_and_canonical_order(self) -> None:
        # POL then GAG genomically, main order POL,GAG -> canonical, 2/3 complete.
        loci = [self._locus({"POL": (100, 200), "GAG": (210, 300)})]
        hits = {
            "L0|POL": [("Gammaretrovirus", 100.0)],
            "L0|GAG": [("Gammaretrovirus", 90.0)],
        }
        rec = tcl._assemble(
            loci, hits, {}, "vTEST", ["POL", "GAG", "ENV"], {}, 0.10, _AXIS
        )[0]
        assert rec["completeness"] == f"{2 / 3:.3f}"
        assert rec["canonical_order"] == "True"
        assert rec["n_main_genes"] == "2"
        assert rec["taxon_call"] == "Gammaretrovirus"
        assert rec["resolved"] == "True"
        assert rec["is_mosaic"] == "False"

    def test_mosaic_flagged_when_genes_disagree(self) -> None:
        loci = [self._locus({"POL": (100, 200), "GAG": (210, 300)})]
        hits = {
            "L0|POL": [("Gammaretrovirus", 100.0)],
            "L0|GAG": [("Betaretrovirus", 100.0)],
        }
        rec = tcl._assemble(loci, hits, {}, "vTEST", ["POL", "GAG"], {}, 0.10, _AXIS)[0]
        assert rec["is_mosaic"] == "True"
        assert "POL:Gammaretrovirus" in rec["mosaic_composition"]
        assert "GAG:Betaretrovirus" in rec["mosaic_composition"]

    def test_diagnostic_gene_calls_by_presence(self) -> None:
        loci = [self._locus({"REX": (100, 200)})]
        rec = tcl._assemble(
            loci, {}, {}, "vTEST", ["POL"], {"REX": "Deltaretrovirus"}, 0.10, _AXIS
        )[0]
        assert rec["taxon_call"] == "Deltaretrovirus"
        assert rec["method"] == "presence"


class TestAxisResolution:
    """ADR-008: 'resolved' = the call landed on a declared axis taxon, at whatever
    rank, replacing the old rank=='genus' test. Off-axis LCA-backoffs are recorded
    (taxon_call + rank) but not marked resolved."""

    def _locus(self, genes: dict[str, tuple[int, int]]) -> dict:
        return {
            "id": "L0",
            "seqname": "chr1",
            "parent": "r",
            "strand": "+",
            "start": min(s for s, _ in genes.values()),
            "end": max(e for _, e in genes.values()),
            "genes": genes,
            "probe_label_set": "",
        }

    def test_axis_member_is_resolved(self) -> None:
        loci = [self._locus({"POL": (100, 200)})]
        hits = {"L0|POL": [("Lentivirus", 100.0)]}
        rec = tcl._assemble(loci, hits, {}, "v", ["POL"], {}, 0.10, {"Lentivirus"})[0]
        assert rec["taxon_call"] == "Lentivirus"
        assert rec["resolved"] == "True"

    def test_backoff_off_axis_is_not_resolved(self) -> None:
        # Near-tie Gamma/Beta -> weighted-LCA backs off to Orthoretrovirinae
        # (a subfamily, NOT an axis member) -> recorded but resolved == False.
        loci = [self._locus({"POL": (100, 200)})]
        hits = {"L0|POL": [("Gammaretrovirus", 600.0), ("Betaretrovirus", 590.0)]}
        rec = tcl._assemble(loci, hits, {}, "v", ["POL"], {}, 0.10, _AXIS)[0]
        assert rec["taxon_call"] == "Orthoretrovirinae"
        assert rec["rank"] == "subfamily"
        assert rec["resolved"] == "False"

    def test_non_genus_axis_member_resolves(self, tmp_path) -> None:
        # A family-rank axis member (Bornaviridae) is a first-class resolved call —
        # the rank-agnostic guarantee. Requires the taxon in the loaded taxonomy.
        tsv = tmp_path / "taxonomy.tsv"
        tsv.write_text(
            "name\tparent\trank\nRiboviria\t\trealm\nBornaviridae\tRiboviria\tfamily\n",
            encoding="utf-8",
        )
        saved_parent, saved_rank = dict(tlca.RETRO_PARENT), dict(tlca.RANK_OF)
        try:
            tlca.load_taxonomy(tsv)
            loci = [self._locus({"POL": (100, 200)})]
            hits = {"L0|POL": [("Bornaviridae", 100.0)]}
            rec = tcl._assemble(
                loci, hits, {}, "v", ["POL"], {}, 0.10, {"Bornaviridae"}
            )[0]
            assert rec["taxon_call"] == "Bornaviridae"
            assert rec["rank"] == "family"
            assert rec["resolved"] == "True"
        finally:
            tlca.RETRO_PARENT, tlca.RANK_OF = saved_parent, saved_rank


class TestGappaParse:
    def test_parses_taxopath_and_confidence(self, tmp_path) -> None:
        tsv = tmp_path / "per_query.tsv"
        tsv.write_text(
            "name\tLWR\taLWR\ttaxopath\n"
            "L0|POL\t0.9\t0.95\tRetroviridae;Orthoretrovirinae;Gammaretrovirus\n",
            encoding="utf-8",
        )
        parsed = tplace._parse_gappa(tsv)
        path, conf = parsed["L0|POL"]
        assert path.endswith("Gammaretrovirus")
        assert abs(conf - 0.95) < 1e-9

    def test_missing_file_returns_empty(self, tmp_path) -> None:
        assert tplace._parse_gappa(tmp_path / "nope.tsv") == {}


class TestPlaceableQuery:
    """A placement query must overlap the reference at >= _MIN_PLACEMENT_SITES
    standard-AA columns; degenerate rows (all-gap, all-X, or a single residue)
    crash epa-ng ('no non-gap sites'), so they are dropped and fall back to LCA."""

    def test_informative_site_count_ignores_gaps_and_ambiguous(self) -> None:
        assert tplace._informative_site_count("--M-K-P--") == 3
        assert tplace._informative_site_count("acdefg") == 6  # lowercase counts
        assert tplace._informative_site_count("---X-*.?-") == 0

    def test_enough_sites_is_placeable(self) -> None:
        assert tplace._is_placeable("A" * tplace._MIN_PLACEMENT_SITES)
        assert tplace._is_placeable(
            "-A-" * tplace._MIN_PLACEMENT_SITES
        )  # gaps interspersed

    def test_degenerate_rows_dropped(self) -> None:
        assert not tplace._is_placeable("---------")  # all gap
        assert not tplace._is_placeable("XXXX")  # translated all-stop -> X
        assert not tplace._is_placeable(
            "-" * 2000 + "F" + "-" * 800
        )  # the L8113 case: 1 site
        assert not tplace._is_placeable(
            "A" * (tplace._MIN_PLACEMENT_SITES - 1)
        )  # just under
        assert not tplace._is_placeable("")


class TestConfidenceTag:
    """confidence_tag (HC/LC) is derived from the locus confidence vs a
    user-adjustable threshold (classification.confidence_min, default 0.5),
    for every method alike."""

    def _locus(self, genes: dict[str, tuple[int, int]]) -> dict:
        return {
            "id": "L0",
            "seqname": "chr1",
            "parent": "r",
            "strand": "+",
            "start": min(s for s, _ in genes.values()),
            "end": max(e for _, e in genes.values()),
            "genes": genes,
            "probe_label_set": "",
        }

    def test_clean_call_is_high_confidence(self) -> None:
        loci = [self._locus({"POL": (100, 200)})]
        hits = {"L0|POL": [("Lentivirus", 100.0)]}  # single taxon -> conf 1.000
        rec = tcl._assemble(loci, hits, {}, "v", ["POL"], {}, 0.10, _AXIS)[0]
        assert rec["confidence"] == "1.000"
        assert rec["confidence_tag"] == "HC"

    def test_unclassified_is_low_confidence(self) -> None:
        loci = [self._locus({"POL": (100, 200)})]
        rec = tcl._assemble(loci, {}, {}, "v", ["POL"], {}, 0.10, _AXIS)[0]  # no hits
        assert rec["taxon_call"] == tlca.UNCLASSIFIED
        assert rec["confidence"] == "0.000"
        assert rec["confidence_tag"] == "LC"

    def test_threshold_is_inclusive_at_boundary(self) -> None:
        # presence call has confidence exactly 1.000; threshold 1.0 -> still HC
        loci = [self._locus({"REX": (100, 200)})]
        rec = tcl._assemble(
            loci,
            {},
            {},
            "v",
            ["POL"],
            {"REX": "Deltaretrovirus"},
            0.10,
            _AXIS,
            confidence_min=1.0,
        )[0]
        assert rec["confidence"] == "1.000"
        assert rec["confidence_tag"] == "HC"

    def test_threshold_is_config_adjustable(self) -> None:
        # a clean 1.000 call tagged LC only under an (extreme) threshold above 1
        loci = [self._locus({"POL": (100, 200)})]
        hits = {"L0|POL": [("Lentivirus", 100.0)]}
        rec = tcl._assemble(
            loci, hits, {}, "v", ["POL"], {}, 0.10, _AXIS, confidence_min=1.5
        )[0]
        assert rec["confidence_tag"] == "LC"


class TestBlastxEvidenceAndSource:
    def _locus(self, genes: dict[str, tuple[int, int]]) -> dict:
        return {
            "id": "L0",
            "seqname": "chr1",
            "parent": "r",
            "strand": "+",
            "start": min(s for s, _ in genes.values()),
            "end": max(e for _, e in genes.values()),
            "genes": genes,
            "probe_label_set": "",
        }

    def test_n_blastx_hits_counts_all_gene_hits(self) -> None:
        loci = [self._locus({"POL": (100, 200), "GAG": (210, 300)})]
        hits = {
            "L0|POL": [("Gammaretrovirus", 100.0), ("Gammaretrovirus", 90.0)],
            "L0|GAG": [("Gammaretrovirus", 80.0)],
        }
        rec = tcl._assemble(loci, hits, {}, "v", ["POL", "GAG"], {}, 0.10, _AXIS)[0]
        assert rec["n_blastx_hits"] == "3"

    def test_no_blastx_hit_is_zero(self) -> None:
        # the candidate-novel-retrovirus signal: valid structure, zero homology
        loci = [self._locus({"POL": (100, 200)})]
        rec = tcl._assemble(loci, {}, {}, "v", ["POL"], {}, 0.10, _AXIS)[0]
        assert rec["n_blastx_hits"] == "0"

    def test_source_defaults_anchored(self) -> None:
        loci = [self._locus({"POL": (100, 200)})]
        rec = tcl._assemble(loci, {}, {}, "v", ["POL"], {}, 0.10, _AXIS)[0]
        assert rec["source"] == "anchored"

    def test_source_is_stamped(self) -> None:
        loci = [self._locus({"POL": (100, 200)})]
        rec = tcl._assemble(
            loci, {}, {}, "v", ["POL"], {}, 0.10, _AXIS, source="fragment"
        )[0]
        assert rec["source"] == "fragment"


class TestGateAndCounts:
    def _rec(self, taxon: str, n_hits: str) -> dict[str, str]:
        return {"taxon_call": taxon, "n_blastx_hits": n_hits}

    def test_gate_drops_only_unclassified(self) -> None:
        recs = [
            self._rec("Lentivirus", "5"),
            self._rec(tlca.UNCLASSIFIED, "0"),
            self._rec("Gammaretrovirus", "3"),
        ]
        kept = tcl.gate_classified(recs)
        assert [r["taxon_call"] for r in kept] == ["Lentivirus", "Gammaretrovirus"]

    def test_anchored_counts(self) -> None:
        recs = [
            self._rec("Lentivirus", "5"),
            self._rec(tlca.UNCLASSIFIED, "0"),
            self._rec(tlca.UNCLASSIFIED, "2"),
        ]
        counts = {
            c["metric"]: c["value"]
            for c in tcl.classification_counts(recs, recs, "anchored")
        }
        assert counts["loci_total"] == 3
        assert counts["loci_classified"] == 1
        assert counts["loci_unclassified"] == 2
        assert counts["loci_no_blastx_hit"] == 1  # only the n_hits == "0" one

    def test_fragment_counts_use_pre_and_post_gate(self) -> None:
        all_recs = [
            self._rec("Lentivirus", "5"),
            self._rec(tlca.UNCLASSIFIED, "0"),
        ]
        kept = tcl.gate_classified(all_recs)
        counts = {
            c["metric"]: c["value"]
            for c in tcl.classification_counts(all_recs, kept, "fragment")
        }
        assert counts["fragments_total"] == 2
        assert counts["fragments_recovered"] == 1


def test_loci_columns_cover_record_keys() -> None:
    """The fixed parquet schema must include every key the assembler emits."""
    loci = [
        {
            "id": "L0",
            "seqname": "chr1",
            "parent": "r",
            "strand": "+",
            "start": 1,
            "end": 9,
            "genes": {"POL": (1, 9)},
            "probe_label_set": "",
        }
    ]
    rec = tcl._assemble(
        loci, {"L0|POL": [("Lentivirus", 50.0)]}, {}, "v", ["POL"], {}, 0.10, _AXIS
    )[0]
    assert set(rec).issubset(set(tcl.LOCI_COLUMNS))
    # erv_class resolves through the loaded/fallback map
    assert "erv_class" in rec
    assert rec["taxon_call"] == "Lentivirus"
    assert rec["resolved"] == "True"
    assert tlca.rank_of("Lentivirus") == "genus"
