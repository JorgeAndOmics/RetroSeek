"""Unit tests for ``workflow/scripts/solo_ltr_integrator.py``.

Covers the four parsers (``parse_valid_ranges``, ``parse_ltr_library_headers``,
``parse_nmtf_pass_list``, ``_parse_coord_column``) and the two label-
propagation paths (consensus-family primary + nearest-ERV fallback) plus
the solo/intact ratio computation.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

# solo_ltr_integrator imports pandas at module load. Skip the whole file
# in environments where pandas isn't installed (the retroseek conda env
# always has it).
pytest.importorskip("pandas")

from solo_ltr_integrator import (
    SoloLTR,
    ValidRange,
    _parse_coord_column,
    _parse_probes_from_gff3_attrs,
    compute_solo_intact_ratio,
    parse_ltr_library_headers,
    parse_nmtf_pass_list,
    parse_valid_ranges,
    propagate_labels,
)


# ---------------------------------------------------------------------
# _parse_probes_from_gff3_attrs
# ---------------------------------------------------------------------
def test_parse_probes_handles_probe_labels_attribute() -> None:
    """``probe_labels=`` is the canonical multi-label encoding."""
    attrs = "ID=erv1;probe_labels=ENV,GAG;source=ranges_analysis"
    assert _parse_probes_from_gff3_attrs(attrs) == ["ENV", "GAG"]


def test_parse_probes_handles_legacy_probe_attribute() -> None:
    """``probe=`` (singular) is accepted as a fallback."""
    attrs = "ID=erv2;probe=POL"
    assert _parse_probes_from_gff3_attrs(attrs) == ["POL"]


def test_parse_probes_returns_empty_when_no_probe_attr() -> None:
    assert _parse_probes_from_gff3_attrs("ID=erv3;other=value") == []


# ---------------------------------------------------------------------
# parse_valid_ranges
# ---------------------------------------------------------------------
def test_parse_valid_ranges_extracts_id_probes_coords(tmp_path: Path) -> None:
    """ID + probe_labels + coords (1→0-indexed) are read off the GFF3 row."""
    gff = tmp_path / "valid.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tRS\tERV\t101\t300\t.\t+\t.\tID=erv1;probe_labels=ENV,GAG\n"
    )
    ranges = parse_valid_ranges(gff)
    assert len(ranges) == 1
    r = ranges[0]
    assert r.chrom == "chr1"
    assert r.start == 100  # 1-indexed 101 → 0-indexed 100
    assert r.end == 300
    assert r.probes == ["ENV", "GAG"]
    assert r.erv_id == "erv1"


# ---------------------------------------------------------------------
# parse_ltr_library_headers
# ---------------------------------------------------------------------
def test_parse_ltr_library_headers_with_source_token(tmp_path: Path) -> None:
    """A header carrying ``source=ID,ID,...`` yields the family→sources map."""
    fa = tmp_path / "lib.fa"
    fa.write_text(
        ">family1#LTR/unknown source=erv1,erv2,erv3\n"
        "ACGTACGTACGT\n"
        ">family2#LTR/Copia members=erv4,erv5\n"
        "TTTTAAAA\n"
    )
    mapping = parse_ltr_library_headers(fa)
    assert mapping == {
        "family1": ["erv1", "erv2", "erv3"],
        "family2": ["erv4", "erv5"],
    }


def test_parse_ltr_library_headers_no_source_returns_empty_list(tmp_path: Path) -> None:
    """A header without ``source=`` or ``members=`` maps to an empty list.

    Triggers the nearest-ERV fallback in propagate_labels.
    """
    fa = tmp_path / "lib_nosource.fa"
    fa.write_text(">family99#LTR/unknown\nACGT\n")
    mapping = parse_ltr_library_headers(fa)
    assert mapping == {"family99": []}


def test_parse_ltr_library_headers_missing_file_returns_empty(tmp_path: Path) -> None:
    """Missing LTRlib.fa is non-fatal — fallback path covers the case."""
    assert parse_ltr_library_headers(tmp_path / "missing.fa") == {}


# ---------------------------------------------------------------------
# _parse_coord_column / parse_nmtf_pass_list
# ---------------------------------------------------------------------
def test_parse_coord_column_basic() -> None:
    assert _parse_coord_column("chr1:100..500") == ("chr1", 100, 500)


def test_parse_coord_column_with_strand_suffix() -> None:
    """``chr1:100..500+`` and ``chr1:100..500-`` both parse the coords."""
    assert _parse_coord_column("chr1:100..500+") == ("chr1", 100, 500)
    assert _parse_coord_column("chr1:100..500-") == ("chr1", 100, 500)


def test_parse_coord_column_rejects_unparseable() -> None:
    assert _parse_coord_column("not a coord") is None


def test_parse_nmtf_pass_list_extracts_coords_and_family(tmp_path: Path) -> None:
    nmtf = tmp_path / "tiny.nmtf.pass.list"
    nmtf.write_text(
        "# header\nchr1:100..500+ family1 95.0\nchr1:1000..1500- family2 88.0\n"
    )
    solos = parse_nmtf_pass_list(nmtf)
    assert len(solos) == 2
    assert (solos[0].chrom, solos[0].start, solos[0].end) == ("chr1", 100, 500)
    assert solos[0].strand == "+"
    assert solos[0].family == "family1"
    assert solos[1].strand == "-"


# ---------------------------------------------------------------------
# propagate_labels — primary path (family) and fallback (nearest_erv)
# ---------------------------------------------------------------------
def test_propagate_labels_uses_family_path_when_resolvable() -> None:
    """When source-ERV IDs map to valid_ranges entries, label_source=family."""
    valid = [
        ValidRange(chrom="chr1", start=99, end=300, probes=["ENV"], erv_id="erv1"),
        ValidRange(chrom="chr1", start=999, end=1200, probes=["GAG"], erv_id="erv2"),
    ]
    solos = [SoloLTR(chrom="chr1", start=5000, end=5100, family="family_A")]
    family_to_sources = {"family_A": ["erv1", "erv2"]}

    propagate_labels(solos, valid, family_to_sources, max_distance=10000)

    assert solos[0].label_source == "family"
    assert sorted(solos[0].probe_labels) == ["ENV", "GAG"]
    assert sorted(solos[0].contributing_ervs) == ["erv1", "erv2"]


def test_propagate_labels_falls_back_to_nearest_erv() -> None:
    """When no source IDs match valid_ranges, take the nearest valid ERV's labels."""
    valid = [
        ValidRange(chrom="chr1", start=4900, end=4950, probes=["POL"], erv_id="ervN"),
    ]
    solos = [SoloLTR(chrom="chr1", start=5000, end=5100, family="family_X")]
    family_to_sources = {"family_X": ["unknown_id"]}

    propagate_labels(solos, valid, family_to_sources, max_distance=200)

    assert solos[0].label_source == "nearest_erv"
    assert solos[0].probe_labels == ["POL"]
    assert solos[0].contributing_ervs == ["ervN"]


def test_propagate_labels_no_label_when_nearest_outside_window() -> None:
    """Distance > max_distance → label_source stays ``none``."""
    valid = [
        ValidRange(chrom="chr1", start=0, end=100, probes=["POL"], erv_id="erv_far"),
    ]
    solos = [SoloLTR(chrom="chr1", start=50000, end=50100, family="x")]
    propagate_labels(solos, valid, {}, max_distance=1000)
    assert solos[0].label_source == "none"
    assert solos[0].probe_labels == []


# ---------------------------------------------------------------------
# compute_solo_intact_ratio
# ---------------------------------------------------------------------
def test_compute_solo_intact_ratio_exclusive_vs_shared_modes() -> None:
    """Single-label rows go to ``exclusive``; multi-label to ``shared``."""
    valid = [
        ValidRange(chrom="chr1", start=0, end=100, probes=["ENV"], erv_id="e1"),
        ValidRange(
            chrom="chr1", start=200, end=300, probes=["ENV", "GAG"], erv_id="e2"
        ),
    ]
    solos = [SoloLTR(chrom="chr1", start=10, end=20, probe_labels=["ENV"])]

    df = compute_solo_intact_ratio(solos, valid, species="Toyus")

    rows = {(r.probe_family, r.label_mode): r for r in df.itertuples()}
    assert rows[("ENV", "exclusive")].intact_count == 1
    assert rows[("ENV", "exclusive")].solo_count == 1
    assert rows[("ENV", "shared")].intact_count == 1
    assert rows[("GAG", "shared")].intact_count == 1


def test_compute_solo_intact_ratio_zero_intact_yields_nan() -> None:
    """When a probe has solos but no intact ERVs, the ratio is NaN."""
    valid: list[ValidRange] = []
    solos = [SoloLTR(chrom="chr1", start=10, end=20, probe_labels=["ENV"])]
    df = compute_solo_intact_ratio(solos, valid, species="Toyus")
    assert len(df) == 1
    row = df.iloc[0]
    assert row["intact_count"] == 0
    assert row["solo_count"] == 1
    assert math.isnan(row["solo_to_intact_ratio"])
