"""Unit tests for the per-locus classifier internals (taxonomy_classify_loci).

These pin the data-derived, gene-agnostic contract: diagnostic genes are
discovered from the reference (never hard-coded), loci are grouped by their
LTR-element ``Parent``, and the assembled record carries the genus-founded
structural metrics (completeness, canonical order, mosaic). A separate test
covers the gappa-output parser used by the placement branch.
"""

from __future__ import annotations

import taxonomy_classify_loci as tcl
import taxonomy_lca as tlca
import taxonomy_placement as tplace


class TestAutoDiagnostic:
    def test_single_genus_gene_is_diagnostic(self) -> None:
        gene_of = {"a1": "REX", "a2": "REX", "a3": "POL"}
        genus_of = {
            "a1": "Deltaretrovirus",
            "a2": "Deltaretrovirus",
            "a3": "Gammaretrovirus",
        }
        diag = tcl.auto_diagnostic(gene_of, genus_of)
        assert diag["REX"] == "Deltaretrovirus"

    def test_multi_genus_gene_not_diagnostic(self) -> None:
        gene_of = {"a1": "POL", "a2": "POL"}
        genus_of = {"a1": "Gammaretrovirus", "a2": "Betaretrovirus"}
        assert "POL" not in tcl.auto_diagnostic(gene_of, genus_of)

    def test_other_catchall_never_diagnostic(self) -> None:
        gene_of = {"a1": "OTHER"}
        genus_of = {"a1": "Lentivirus"}
        assert "OTHER" not in tcl.auto_diagnostic(gene_of, genus_of)


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
        rec = tcl._assemble(loci, hits, {}, "vTEST", ["POL", "GAG", "ENV"], {}, 0.10)[0]
        assert rec["completeness"] == f"{2 / 3:.3f}"
        assert rec["canonical_order"] == "True"
        assert rec["n_main_genes"] == "2"
        assert rec["genus_call"] == "Gammaretrovirus"
        assert rec["is_mosaic"] == "False"

    def test_mosaic_flagged_when_genes_disagree(self) -> None:
        loci = [self._locus({"POL": (100, 200), "GAG": (210, 300)})]
        hits = {
            "L0|POL": [("Gammaretrovirus", 100.0)],
            "L0|GAG": [("Betaretrovirus", 100.0)],
        }
        rec = tcl._assemble(loci, hits, {}, "vTEST", ["POL", "GAG"], {}, 0.10)[0]
        assert rec["is_mosaic"] == "True"
        assert "POL:Gammaretrovirus" in rec["mosaic_composition"]
        assert "GAG:Betaretrovirus" in rec["mosaic_composition"]

    def test_diagnostic_gene_calls_by_presence(self) -> None:
        loci = [self._locus({"REX": (100, 200)})]
        rec = tcl._assemble(
            loci, {}, {}, "vTEST", ["POL"], {"REX": "Deltaretrovirus"}, 0.10
        )[0]
        assert rec["genus_call"] == "Deltaretrovirus"
        assert rec["method"] == "presence"


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


class TestInformativeResidues:
    """A placement query must carry >=1 standard amino acid; an all-gap or
    all-X (degenerate) row crashes epa-ng ('no non-gap sites'), so it is dropped."""

    def test_standard_residues_are_informative(self) -> None:
        assert tplace._has_informative_residues("--M-K-P--")
        assert tplace._has_informative_residues("acdefg")  # lowercase counts

    def test_all_gap_or_all_X_is_dropped(self) -> None:
        assert not tplace._has_informative_residues("---------")
        assert not tplace._has_informative_residues("XXXX")  # translated all-stop -> X
        assert not tplace._has_informative_residues("-X-.X*?-")  # gaps + ambiguous only
        assert not tplace._has_informative_residues("")


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
        hits = {"L0|POL": [("Lentivirus", 100.0)]}  # single genus -> conf 1.000
        rec = tcl._assemble(loci, hits, {}, "v", ["POL"], {}, 0.10)[0]
        assert rec["confidence"] == "1.000"
        assert rec["confidence_tag"] == "HC"

    def test_unclassified_is_low_confidence(self) -> None:
        loci = [self._locus({"POL": (100, 200)})]
        rec = tcl._assemble(loci, {}, {}, "v", ["POL"], {}, 0.10)[0]  # no hits
        assert rec["genus_call"] == tlca.UNCLASSIFIED
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
            confidence_min=1.0,
        )[0]
        assert rec["confidence"] == "1.000"
        assert rec["confidence_tag"] == "HC"

    def test_threshold_is_config_adjustable(self) -> None:
        # a clean 1.000 call tagged LC only under an (extreme) threshold above 1
        loci = [self._locus({"POL": (100, 200)})]
        hits = {"L0|POL": [("Lentivirus", 100.0)]}
        rec = tcl._assemble(loci, hits, {}, "v", ["POL"], {}, 0.10, confidence_min=1.5)[
            0
        ]
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
        rec = tcl._assemble(loci, hits, {}, "v", ["POL", "GAG"], {}, 0.10)[0]
        assert rec["n_blastx_hits"] == "3"

    def test_no_blastx_hit_is_zero(self) -> None:
        # the candidate-novel-retrovirus signal: valid structure, zero homology
        loci = [self._locus({"POL": (100, 200)})]
        rec = tcl._assemble(loci, {}, {}, "v", ["POL"], {}, 0.10)[0]
        assert rec["n_blastx_hits"] == "0"

    def test_source_defaults_anchored(self) -> None:
        loci = [self._locus({"POL": (100, 200)})]
        rec = tcl._assemble(loci, {}, {}, "v", ["POL"], {}, 0.10)[0]
        assert rec["source"] == "anchored"

    def test_source_is_stamped(self) -> None:
        loci = [self._locus({"POL": (100, 200)})]
        rec = tcl._assemble(loci, {}, {}, "v", ["POL"], {}, 0.10, source="fragment")[0]
        assert rec["source"] == "fragment"


class TestGateAndCounts:
    def _rec(self, genus: str, n_hits: str) -> dict[str, str]:
        return {"genus_call": genus, "n_blastx_hits": n_hits}

    def test_gate_drops_only_unclassified(self) -> None:
        recs = [
            self._rec("Lentivirus", "5"),
            self._rec(tlca.UNCLASSIFIED, "0"),
            self._rec("Gammaretrovirus", "3"),
        ]
        kept = tcl.gate_classified(recs)
        assert [r["genus_call"] for r in kept] == ["Lentivirus", "Gammaretrovirus"]

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
        loci, {"L0|POL": [("Lentivirus", 50.0)]}, {}, "v", ["POL"], {}, 0.10
    )[0]
    assert set(rec).issubset(set(tcl.LOCI_COLUMNS))
    # erv_class resolves through the loaded/fallback map
    assert "erv_class" in rec
    assert rec["genus_call"] == "Lentivirus"
    assert tlca.rank_of("Lentivirus") == "genus"
