"""Unit tests for the reference builder's axis resolution (no network).

Pins the ADR-008 hybrid, rank-agnostic axis contract: an explicit ``--taxa`` list
wins; otherwise the axis is derived from the distinct probeset ``Label`` values (any
rank); otherwise the retroviral-genus default. Order-preserving and de-duplicated.
"""

from __future__ import annotations

import taxonomy_reference_builder as trb


class TestResolveAxis:
    def test_explicit_taxa_win(self, tmp_path) -> None:
        csv = tmp_path / "probes.csv"
        csv.write_text("Label,Probe\nLentivirus,POL\n", encoding="utf-8")
        # explicit taxa override the probe CSV entirely
        assert trb.resolve_axis(["Gammaretrovirus", "Bornaviridae"], csv) == [
            "Gammaretrovirus",
            "Bornaviridae",
        ]

    def test_derives_from_probe_label_when_no_taxa(self, tmp_path) -> None:
        csv = tmp_path / "probes.csv"
        # mixed rank (genus + family) + a duplicate -> de-duplicated, order-preserved
        csv.write_text(
            "Label,Abbreviation,Name,Probe,Accession\n"
            "Lentivirus,HIV,x,POL,A1\n"
            "Bornaviridae,BDV,y,N,A2\n"
            "Lentivirus,SIV,z,GAG,A3\n",
            encoding="utf-8",
        )
        assert trb.resolve_axis([], csv) == ["Lentivirus", "Bornaviridae"]

    def test_falls_back_to_default_taxa(self) -> None:
        assert trb.resolve_axis([], None) == list(trb.DEFAULT_TAXA)

    def test_blank_entries_are_dropped(self, tmp_path) -> None:
        assert trb.resolve_axis(["Lentivirus", "", "  "], None) == ["Lentivirus"]
