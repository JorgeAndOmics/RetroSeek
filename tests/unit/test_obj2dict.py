"""Tests for obj2dict.extract_attributes_from_object, the row behind the BLAST table.

The parquet's column order and its None-filling for a missing GenBank,
alignment or HSP record are what the R stages read, so both are pinned here.
"""

from __future__ import annotations

from types import SimpleNamespace

import obj2dict

_BASE = [
    "label", "virus", "abbreviation", "species", "probe", "accession",
    "identifier", "strand", "species_name",
]  # fmt: skip
_GENBANK = [
    "genbank_id", "genbank_name", "genbank_description", "genbank_dbxrefs",
    "genbank_annotations", "genbank_seq",
]  # fmt: skip
_ALIGNMENT = [
    "alignment_title", "alignment_length", "alignment_accession",
    "alignment_hit_id", "alignment_hit_def",
]  # fmt: skip
_HSP = [
    "hsp_bits", "hsp_score", "hsp_evalue", "hsp_query", "hsp_sbjct",
    "hsp_query_start", "hsp_query_end", "hsp_sbjct_start", "hsp_sbjct_end",
    "hsp_identity", "hsp_align_length", "hsp_gaps", "hsp_positives",
    "hsp_strand", "hsp_frame",
]  # fmt: skip


def _probe(**records: object) -> SimpleNamespace:
    base = {
        "label": "ALV", "virus": "Avian leukosis virus", "abbreviation": "ALV",
        "species": "Toyus_toyus", "probe": "POL",
        "accession": "CM1.1 Toyus toyus chromosome 1", "identifier": "abc123",
        "strand": "+", "genbank": None, "alignment": None, "HSP": None,
    }  # fmt: skip
    return SimpleNamespace(**(base | records))


def test_columns_come_in_a_fixed_order() -> None:
    row = obj2dict.extract_attributes_from_object(_probe())
    assert list(row) == _BASE + _GENBANK + _ALIGNMENT + _HSP


def test_missing_records_fill_their_columns_with_none() -> None:
    row = obj2dict.extract_attributes_from_object(_probe())
    assert all(row[c] is None for c in _GENBANK + _ALIGNMENT + _HSP)


def test_accession_keeps_only_the_first_token() -> None:
    row = obj2dict.extract_attributes_from_object(_probe())
    assert row["accession"] == "CM1.1"


def test_present_records_are_read() -> None:
    gb = SimpleNamespace(
        id="A1", name="n", description="d", dbxrefs=["x"], annotations={"k": 1},
        seq="MKV",
    )  # fmt: skip
    hsp = SimpleNamespace(
        bits=55.1, score=120, expect=1e-10, query="MK", sbjct="MK", query_start=1,
        query_end=2, sbjct_start=100, sbjct_end=105, identities=2, align_length=2,
        gaps=0, positives=2, strand=(None, None), frame=(0, 2),
    )  # fmt: skip
    row = obj2dict.extract_attributes_from_object(_probe(genbank=gb, HSP=hsp))
    assert row["genbank_annotations"] == "{'k': 1}"
    assert row["genbank_seq"] == "MKV"
    assert (row["hsp_evalue"], row["hsp_identity"], row["hsp_frame"]) == (
        1e-10,
        2,
        (0, 2),
    )
    assert row["alignment_title"] is None
