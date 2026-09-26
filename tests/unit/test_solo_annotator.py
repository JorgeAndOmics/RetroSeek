"""Unit tests for ``workflow/scripts/solo_ltr/solo_annotator.py``.

The integrator annotates LTR_retriever's solo LTRs with RetroSeek's own
taxonomy. Two contracts are pinned here.

**Where solos come from.** ``solo_finder.pl`` writes
``chrom, start, end, locus, library_id, coverage``, driven off the whole-genome
RepeatMasker table. The retired implementation read ``nmtf.pass.list``, which
holds *intact* non-TGCA elements, and would have reported them as solos.

**How a solo inherits a taxon.** LTR_retriever names each library sequence
after the genomic span of the intact element that seeded it
(``>{chr}:{start}..{end}#LTR/{fam}``, annotate_lib.pl), so the library ID is a
coordinate, not an opaque ``family1``. Overlapping that span against the
classified loci table inherits ``taxon_call`` from the element the solo's
sequence actually came from. Overlap rather than exact equality, because
LTR_retriever adjusts element boundaries, and a locus sits inside its element
rather than sharing its edges.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from solo_annotator import (
    annotate_solos,
    compute_solo_intact_ratio,
    parse_library_locus,
    parse_loci_csv,
    parse_solo_list,
    write_solo_ltr_gff3,
    write_solo_table,
)

LOCI_HEADER = (
    "id,seqname,start,end,strand,parent,genes_present,n_main_genes,completeness,"
    "canonical_order,structure_class,domain_tier,oversized,taxon_call,rank,segment,"
    "segment_rank,resolved,confidence,confidence_tag,n_blastx_hits,method,per_gene,"
    "is_mosaic,mosaic_composition,erv_class,probe_label_set,ref_version,source"
)


def _locus_row(
    locus_id: str,
    seqname: str,
    start: int,
    end: int,
    taxon: str,
    segment: str,
    source: str = "ltr-flanked",
    erv_class: str = "Class I",
    confidence: str = "1.000",
) -> str:
    return (
        f'{locus_id},{seqname},{start},{end},+,LTR_retrotransposon1,"POL",1,0.333,'
        f"True,partial,non_domain,False,{taxon},genus,{segment},genus,True,"
        f"{confidence},HC,10,lca,,False,,{erv_class},POL,1,{source}"
    )


@pytest.fixture
def loci_csv(tmp_path: Path) -> Path:
    """Two ltr-flanked loci plus one orphan, which must not act as a donor."""
    path = tmp_path / "Toyus.loci.csv"
    path.write_text(
        LOCI_HEADER
        + "\n"
        + _locus_row("L0", "chr1", 1200, 1400, "Gammaretrovirus", "Gammaretrovirus")
        + "\n"
        + _locus_row("L1", "chr1", 90000, 90500, "Betaretrovirus", "Betaretrovirus")
        + "\n"
        + _locus_row(
            "L2",
            "chr1",
            5000,
            5100,
            "Alpharetrovirus",
            "Alpharetrovirus",
            source="orphan",
        )
        + "\n"
    )
    return path


def _write_solo_list(path: Path, rows: list[tuple[str, int, int, str, float]]) -> Path:
    path.write_text(
        "".join(
            f"{chrom}\t{start}\t{end}\t{chrom}:{start}..{end}\t{library}\t{cov}\n"
            for chrom, start, end, library, cov in rows
        )
    )
    return path


# ---------------------------------------------------------------------
# parse_library_locus
# ---------------------------------------------------------------------
def test_parse_library_locus_reads_chrom_and_span() -> None:
    """``annotate_lib.pl`` writes >{chr}:{start}..{end}#LTR/{fam}."""
    assert parse_library_locus("chr1:1000..5000#LTR/Gypsy") == ("chr1", 1000, 5000)


def test_parse_library_locus_tolerates_a_region_suffix_and_no_family() -> None:
    """Internal-region entries add ``_INT``; RepeatMasker may drop the #class."""
    assert parse_library_locus("chr1:1000..5000_INT#LTR/unknown") == (
        "chr1",
        1000,
        5000,
    )
    assert parse_library_locus("chr1:1000..5000") == ("chr1", 1000, 5000)


def test_parse_library_locus_normalises_reversed_minus_strand_spans() -> None:
    """LTR_retriever writes minus-strand elements start > end.

    Observed on the Desmodus pilot as `CM040301.1:40850808..40843686_LTR`, which
    accounted for 47% of all solos. A reversed interval does not error - the
    shared-base arithmetic just returns zero - so every one of them silently lost
    its taxon and fell through to `label_source=none`.
    """
    assert parse_library_locus("chr1:5000..1000#LTR/Gypsy") == ("chr1", 1000, 5000)


def test_parse_library_locus_returns_none_when_unparseable() -> None:
    """A non-coordinate library name is not an error, just an unusable donor."""
    assert parse_library_locus("family1#LTR/Copia") is None


def test_solo_inherits_taxon_through_a_reversed_library_span(
    tmp_path: Path, loci_csv: Path
) -> None:
    """End-to-end guard: a minus-strand library entry must still find its locus."""
    solo_list = _write_solo_list(
        tmp_path / "s.txt", [("chr1", 7000, 7300, "chr1:5000..1000#LTR/Gypsy", 0.95)]
    )
    solos = parse_solo_list(solo_list)
    annotate_solos(solos, parse_loci_csv(loci_csv), max_distance=10000)

    assert solos[0].taxon_call == "Gammaretrovirus"
    assert solos[0].label_source == "library"


# ---------------------------------------------------------------------
# parse_solo_list / parse_loci_csv
# ---------------------------------------------------------------------
def test_parse_solo_list_reads_six_column_rows(tmp_path: Path) -> None:
    path = _write_solo_list(
        tmp_path / "s.txt", [("chr1", 7000, 7300, "chr1:1000..5000#LTR/Gypsy", 0.95)]
    )
    solos = parse_solo_list(path)
    assert len(solos) == 1
    assert (solos[0].chrom, solos[0].start, solos[0].end) == ("chr1", 7000, 7300)
    assert solos[0].coverage == pytest.approx(0.95)


def test_parse_solo_list_of_a_missing_or_empty_file_is_empty(tmp_path: Path) -> None:
    """No solos is a legitimate biological result, not a failure."""
    assert parse_solo_list(tmp_path / "absent.txt") == []
    empty = tmp_path / "empty.txt"
    empty.write_text("")
    assert parse_solo_list(empty) == []


def test_parse_loci_csv_keeps_only_ltr_flanked_donors(loci_csv: Path) -> None:
    """Orphans have no LTR element, so they cannot have seeded a library entry."""
    loci = parse_loci_csv(loci_csv)
    assert sorted(locus.id for locus in loci) == ["L0", "L1"]
    assert loci[0].taxon_call == "Gammaretrovirus"


# ---------------------------------------------------------------------
# annotate_solos
# ---------------------------------------------------------------------
def test_solo_inherits_taxon_from_the_overlapping_library_element(
    tmp_path: Path, loci_csv: Path
) -> None:
    """The primary path: library span overlaps L0, so the solo becomes L0's taxon.

    Note the solo itself sits at chr1:7000-7300, far from L0 at 1200-1400. The
    inheritance follows sequence homology (which library entry matched), not
    proximity - that is the whole point of using the library ID.
    """
    solo_list = _write_solo_list(
        tmp_path / "s.txt", [("chr1", 7000, 7300, "chr1:1000..5000#LTR/Gypsy", 0.95)]
    )
    solos = parse_solo_list(solo_list)
    annotate_solos(solos, parse_loci_csv(loci_csv), max_distance=10000)

    assert solos[0].taxon_call == "Gammaretrovirus"
    assert solos[0].segment == "Gammaretrovirus"
    assert solos[0].label_source == "library"
    assert solos[0].source_loci == ["L0"]


def test_library_element_spanning_two_loci_picks_the_largest_overlap(
    tmp_path: Path,
) -> None:
    """A wide library element can cover more than one locus; the call must be
    deterministic rather than first-wins."""
    loci = tmp_path / "multi.loci.csv"
    loci.write_text(
        LOCI_HEADER
        + "\n"
        + _locus_row("L0", "chr1", 1000, 1050, "Gammaretrovirus", "Gammaretrovirus")
        + "\n"
        + _locus_row("L1", "chr1", 2000, 2900, "Betaretrovirus", "Betaretrovirus")
        + "\n"
    )
    solo_list = _write_solo_list(
        tmp_path / "s.txt", [("chr1", 9000, 9300, "chr1:900..3000#LTR/Gypsy", 0.9)]
    )
    solos = parse_solo_list(solo_list)
    annotate_solos(solos, parse_loci_csv(loci), max_distance=10000)

    assert solos[0].taxon_call == "Betaretrovirus"  # 900 bp beats 51 bp
    assert solos[0].source_loci == ["L1", "L0"]  # all contributors, best first


def test_solo_falls_back_to_the_nearest_locus_when_the_library_id_is_opaque(
    tmp_path: Path, loci_csv: Path
) -> None:
    """Older LTR_retriever names, or RepeatMasker renames, may lose coordinates.

    The fallback is proximity, which is a weaker signal than homology, so it is
    recorded distinctly in ``label_source`` for downstream filtering.
    """
    solo_list = _write_solo_list(
        tmp_path / "s.txt", [("chr1", 1500, 1800, "family7#LTR/Copia", 0.9)]
    )
    solos = parse_solo_list(solo_list)
    annotate_solos(solos, parse_loci_csv(loci_csv), max_distance=10000)

    assert solos[0].taxon_call == "Gammaretrovirus"  # L0 at 1200-1400 is nearest
    assert solos[0].label_source == "nearest_locus"


def test_fallback_respects_the_distance_ceiling(tmp_path: Path, loci_csv: Path) -> None:
    """Beyond the window the solo stays unlabelled rather than guessing."""
    solo_list = _write_solo_list(
        tmp_path / "s.txt", [("chr1", 500000, 500300, "family7#LTR/Copia", 0.9)]
    )
    solos = parse_solo_list(solo_list)
    annotate_solos(solos, parse_loci_csv(loci_csv), max_distance=10000)

    assert solos[0].taxon_call == ""
    assert solos[0].label_source == "none"


def test_library_path_wins_over_a_closer_neighbour(
    tmp_path: Path, loci_csv: Path
) -> None:
    """Homology beats proximity: a solo sitting on top of L1 but whose library
    entry came from L0 inherits L0."""
    solo_list = _write_solo_list(
        tmp_path / "s.txt", [("chr1", 90100, 90300, "chr1:1000..5000#LTR/Gypsy", 0.95)]
    )
    solos = parse_solo_list(solo_list)
    annotate_solos(solos, parse_loci_csv(loci_csv), max_distance=10000)

    assert solos[0].taxon_call == "Gammaretrovirus"
    assert solos[0].label_source == "library"


# ---------------------------------------------------------------------
# compute_solo_intact_ratio
# ---------------------------------------------------------------------
def test_ratio_groups_by_segment_and_divides_by_intact_loci(
    tmp_path: Path, loci_csv: Path
) -> None:
    """solo / intact per group, where intact is the ltr-flanked locus count.

    Two Gammaretrovirus solos against one Gammaretrovirus locus gives 2.0.
    """
    solo_list = _write_solo_list(
        tmp_path / "s.txt",
        [
            ("chr1", 7000, 7300, "chr1:1000..5000#LTR/Gypsy", 0.95),
            ("chr1", 8000, 8300, "chr1:1000..5000#LTR/Gypsy", 0.92),
        ],
    )
    loci = parse_loci_csv(loci_csv)
    solos = parse_solo_list(solo_list)
    annotate_solos(solos, loci, max_distance=10000)

    frame = compute_solo_intact_ratio(solos, loci, species="Toyus", group_by="segment")
    gamma = frame[frame["group"] == "Gammaretrovirus"].iloc[0]
    assert gamma["solo_count"] == 2
    assert gamma["intact_count"] == 1
    assert gamma["solo_to_intact_ratio"] == pytest.approx(2.0)
    assert gamma["species"] == "Toyus"


def test_ratio_reports_groups_with_intact_loci_but_no_solos(
    tmp_path: Path, loci_csv: Path
) -> None:
    """A zero-solo lineage is a real finding and must appear, not vanish."""
    solo_list = _write_solo_list(
        tmp_path / "s.txt", [("chr1", 7000, 7300, "chr1:1000..5000#LTR/Gypsy", 0.95)]
    )
    loci = parse_loci_csv(loci_csv)
    solos = parse_solo_list(solo_list)
    annotate_solos(solos, loci, max_distance=10000)

    frame = compute_solo_intact_ratio(solos, loci, species="Toyus", group_by="segment")
    beta = frame[frame["group"] == "Betaretrovirus"].iloc[0]
    assert beta["solo_count"] == 0
    assert beta["intact_count"] == 1
    assert beta["solo_to_intact_ratio"] == pytest.approx(0.0)


def test_ratio_group_by_none_pools_every_locus(tmp_path: Path, loci_csv: Path) -> None:
    solo_list = _write_solo_list(
        tmp_path / "s.txt", [("chr1", 7000, 7300, "chr1:1000..5000#LTR/Gypsy", 0.95)]
    )
    loci = parse_loci_csv(loci_csv)
    solos = parse_solo_list(solo_list)
    annotate_solos(solos, loci, max_distance=10000)

    frame = compute_solo_intact_ratio(solos, loci, species="Toyus", group_by="none")
    assert len(frame) == 1
    assert frame.iloc[0]["intact_count"] == 2


def test_ratio_rejects_an_unknown_group_by() -> None:
    with pytest.raises(ValueError, match="group_by"):
        compute_solo_intact_ratio([], [], species="Toyus", group_by="probe_family")


# ---------------------------------------------------------------------
# writers
# ---------------------------------------------------------------------
def test_gff3_carries_taxonomy_attributes(tmp_path: Path, loci_csv: Path) -> None:
    solo_list = _write_solo_list(
        tmp_path / "s.txt", [("chr1", 7000, 7300, "chr1:1000..5000#LTR/Gypsy", 0.95)]
    )
    solos = parse_solo_list(solo_list)
    annotate_solos(solos, parse_loci_csv(loci_csv), max_distance=10000)

    out = tmp_path / "solo.gff3"
    write_solo_ltr_gff3(solos, out, genome="Toyus")
    text = out.read_text()
    assert text.startswith("##gff-version 3")
    assert "solo_LTR" in text
    assert "taxon_call=Gammaretrovirus" in text
    assert "label_source=library" in text
    # The track and the table must share a key, else a feature cannot be
    # looked up in the catalog.
    assert "ID=S0" in text


def test_empty_solo_set_still_writes_a_valid_gff3(tmp_path: Path) -> None:
    """Snakemake needs the file to exist even when a genome has no solos."""
    out = tmp_path / "empty.gff3"
    write_solo_ltr_gff3([], out, genome="Toyus")
    assert out.read_text().startswith("##gff-version 3")


def test_solo_table_uses_catalog_column_names(tmp_path: Path, loci_csv: Path) -> None:
    """The table feeds catalog.csv as the third tier, so it must speak its schema."""
    solo_list = _write_solo_list(
        tmp_path / "s.txt", [("chr1", 7000, 7300, "chr1:1000..5000#LTR/Gypsy", 0.95)]
    )
    solos = parse_solo_list(solo_list)
    annotate_solos(solos, parse_loci_csv(loci_csv), max_distance=10000)

    csv_path = tmp_path / "solo.csv"
    write_solo_table(solos, csv_path, tmp_path / "solo.parquet", species="Toyus rex")
    header = csv_path.read_text().splitlines()[0].split(",")
    for column in (
        "species",
        "source",
        "seqname",
        "start",
        "end",
        "taxon_call",
        "segment",
        "structure_class",
        "id",
    ):
        assert column in header
    body = csv_path.read_text().splitlines()[1]
    assert "solo-ltr" in body
    assert "solo_ltr" in body  # structure_class


# ---------------------------------------------------------------------
# Tie-breaks and malformed rows (pinned in the quality campaign, 2026-09-26)
# ---------------------------------------------------------------------
def test_nearest_locus_tie_goes_to_the_locus_listed_first(tmp_path: Path) -> None:
    loci = tmp_path / "tie.loci.csv"
    loci.write_text(
        LOCI_HEADER
        + "\n"
        + _locus_row("L0", "chr1", 1000, 1099, "Gammaretrovirus", "Gammaretrovirus")
        + "\n"
        + _locus_row("L1", "chr1", 1300, 1399, "Betaretrovirus", "Betaretrovirus")
        + "\n"
    )
    solo_list = _write_solo_list(
        tmp_path / "s.txt", [("chr1", 1150, 1249, "family7#LTR/Copia", 0.9)]
    )
    solos = parse_solo_list(solo_list)
    annotate_solos(solos, parse_loci_csv(loci), max_distance=10000)
    assert solos[0].source_loci == ["L0"]  # both 51 bp away; L0 comes first


def test_equal_library_overlaps_are_ordered_by_locus_id(tmp_path: Path) -> None:
    loci = tmp_path / "tie.loci.csv"
    loci.write_text(
        LOCI_HEADER
        + "\n"
        + _locus_row("L9", "chr1", 1000, 1099, "Gammaretrovirus", "Gammaretrovirus")
        + "\n"
        + _locus_row("L1", "chr1", 2000, 2099, "Betaretrovirus", "Betaretrovirus")
        + "\n"
    )
    solo_list = _write_solo_list(
        tmp_path / "s.txt", [("chr1", 9000, 9300, "chr1:900..3000#LTR/Gypsy", 0.9)]
    )
    solos = parse_solo_list(solo_list)
    annotate_solos(solos, parse_loci_csv(loci), max_distance=10000)
    assert solos[0].source_loci == ["L1", "L9"]  # 100 bp each: id decides
    assert solos[0].taxon_call == "Betaretrovirus"


def test_parse_solo_list_skips_comments_blank_and_unparseable_rows(
    tmp_path: Path,
) -> None:
    path = tmp_path / "s.txt"
    path.write_text(
        "# header\n"
        "\n"
        "\t\t\t\t\t\n"
        "chr1\tx\t20\tchr1:1..20\tlib\t0.5\n"
        "chr1\t1\t20\n"
        "chr1\t1\t20\tchr1:1..20\tlib\t0.5\n"
    )
    solos = parse_solo_list(path)
    assert [(s.chrom, s.start, s.end, s.coverage) for s in solos] == [
        ("chr1", 1, 20, 0.5)
    ]


def test_unparseable_solo_rows_are_reported(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A changed column layout must not pass as 'zero solos' in silence."""
    path = tmp_path / "s.txt"
    path.write_text(
        "chr1\tx\t20\tchr1:1..20\tlib\t0.5\n"
        "chr1\t1\t20\tchr1:1..20\tlib\tnot-a-number\n"
        "chr1\t1\t20\tchr1:1..20\tlib\t0.5\n"
    )
    with caplog.at_level("WARNING"):
        solos = parse_solo_list(path)
    assert len(solos) == 1
    assert "2 of 3 solo rows" in caplog.text


def test_clean_solo_lists_log_no_warning(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    path = _write_solo_list(tmp_path / "s.txt", [("chr1", 1, 20, "lib", 0.5)])
    with caplog.at_level("WARNING"):
        parse_solo_list(path)
    assert caplog.text == ""
