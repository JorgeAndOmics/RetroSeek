"""Tests for the native solo-LTR detector.

The detector is a funnel followed by a subtraction, and both halves can be wrong
quietly. A threshold applied with the wrong comparison changes the solo count by
thousands without erroring; a classification priority applied in the wrong order
turns intact elements into solos. Each test below pins one of those.

Fixture style follows tests/unit/test_scan_domains.py: an inline tabular constant
in the real output format, written to tmp_path by a fixture.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest
import solo_finder
from solo_intervals import IntervalIndex

# blastn -outfmt "6 qseqid sseqid pident length qlen sstart send evalue bitscore".
# One row per case, named in the comment so a failure is readable.
_HITS = """\
chr1|LTR_retrotransposon1|L\tchr1\t99.0\t400\t400\t50000\t50399\t1e-90\t700
chr1|LTR_retrotransposon1|L\tchr1\t99.0\t400\t400\t60000\t60399\t1e-90\t700
chr1|LTR_retrotransposon1|L\tchr1\t70.0\t400\t400\t70000\t70399\t1e-20\t200
chr1|LTR_retrotransposon1|L\tchr1\t99.0\t50\t400\t80000\t80049\t1e-10\t90
chr1|LTR_retrotransposon1|L\tchr1\t99.0\t700\t400\t90000\t90699\t1e-90\t700
"""

# GFF3 of intact elements. Element 1 sits at 49,000..61,000, so it covers the
# first two hits above.
_ELEMENTS = """\
##gff-version 3
chr1\tLTRdigest\tLTR_retrotransposon\t49000\t61000\t.\t+\t.\tID=LTR_retrotransposon1
"""

# An orphan locus 5 kb from the hit at 90,000.
_ORPHANS = """\
##gff-version 3
chr1\trtracklayer\tproviral_sequence\t95000\t95500\t.\t+\t.\tprobe=POL
"""


@pytest.fixture
def thresholds() -> solo_finder.Thresholds:
    """The ADR-017 operating point."""
    return solo_finder.Thresholds(
        min_identity=95.0,
        min_coverage=0.8,
        max_coverage=1.2,
        min_alignment_length=80,
        min_hit_length=300,
        merge_gap=0,
        orphan_pad=10000,
    )


@pytest.fixture
def hits(tmp_path: Path) -> Path:
    p = tmp_path / "hits.tsv.gz"
    with gzip.open(p, "wt") as handle:
        handle.write(_HITS)
    return p


@pytest.fixture
def elements(tmp_path: Path) -> Path:
    p = tmp_path / "elements.gff3"
    p.write_text(_ELEMENTS)
    return p


@pytest.fixture
def orphans(tmp_path: Path) -> Path:
    p = tmp_path / "orphans.gff3"
    p.write_text(_ORPHANS)
    return p


# ---- the acceptance criteria, one test each ----


def test_a_hit_below_the_identity_floor_is_rejected(
    thresholds: solo_finder.Thresholds,
) -> None:
    hit = solo_finder.Hit("bait", "chr1", 70.0, 400, 400, 1, 400)
    assert solo_finder.rejection_reason(hit, thresholds) == "identity"


def test_a_hit_below_the_absolute_alignment_minimum_is_rejected(
    thresholds: solo_finder.Thresholds,
) -> None:
    hit = solo_finder.Hit("bait", "chr1", 99.0, 50, 400, 1, 50)
    # 50 bp fails the 80 bp absolute floor before coverage is even considered.
    assert solo_finder.rejection_reason(hit, thresholds) == "alignment_length"


def test_a_hit_covering_too_little_of_the_bait_is_rejected(
    thresholds: solo_finder.Thresholds,
) -> None:
    hit = solo_finder.Hit("bait", "chr1", 99.0, 200, 400, 1, 200)  # coverage 0.5
    assert solo_finder.rejection_reason(hit, thresholds) == "coverage"


def test_a_hit_running_far_past_the_bait_is_rejected(
    thresholds: solo_finder.Thresholds,
) -> None:
    hit = solo_finder.Hit("bait", "chr1", 99.0, 700, 400, 1, 700)  # coverage 1.75
    assert solo_finder.rejection_reason(hit, thresholds) == "coverage"


def test_a_short_but_proportionate_hit_is_rejected_by_hit_length(
    thresholds: solo_finder.Thresholds,
) -> None:
    """A 100 bp hit covering all of a 100 bp bait passes coverage and fails length.

    This is the criterion that took the prototype from 492:1 to 27.6:1, so it is
    tested separately from coverage even though both are length-ish rules.
    """
    hit = solo_finder.Hit("bait", "chr1", 99.0, 100, 100, 1, 100)
    assert solo_finder.rejection_reason(hit, thresholds) == "hit_length"


def test_a_hit_at_every_boundary_exactly_is_accepted(
    thresholds: solo_finder.Thresholds,
) -> None:
    """Every threshold is inclusive, so the boundary value passes."""
    hit = solo_finder.Hit(
        "bait", "chr1", 95.0, 300, 375, 1, 300
    )  # coverage exactly 0.8
    assert solo_finder.rejection_reason(hit, thresholds) is None


# ---- merging ----


def test_overlapping_hits_become_one_candidate(
    thresholds: solo_finder.Thresholds,
) -> None:
    accepted = [
        solo_finder.Hit("b1", "chr1", 99.0, 400, 400, 1000, 1399),
        solo_finder.Hit("b2", "chr1", 97.0, 400, 400, 1200, 1599),
    ]
    candidates = solo_finder.merge_candidates(accepted, gap=0)
    assert len(candidates) == 1
    assert (candidates[0].start, candidates[0].end) == (1000, 1599)
    assert candidates[0].n_hits == 2


def test_the_candidate_keeps_the_best_identity_and_its_bait(
    thresholds: solo_finder.Thresholds,
) -> None:
    """The representative bait is the best-matching hit, not the first seen.

    It decides which element the solo inherits taxonomy from.
    """
    accepted = [
        solo_finder.Hit("weak", "chr1", 96.0, 400, 400, 1000, 1399),
        solo_finder.Hit("strong", "chr1", 99.5, 400, 400, 1100, 1499),
    ]
    candidate = solo_finder.merge_candidates(accepted, gap=0)[0]
    assert candidate.best_identity == 99.5
    assert candidate.bait == "strong"


def test_an_equal_identity_keeps_the_bait_seen_first() -> None:
    accepted = [
        solo_finder.Hit("first", "chr1", 99.0, 400, 400, 1000, 1399),
        solo_finder.Hit("second", "chr1", 99.0, 400, 400, 1100, 1499),
    ]
    candidate = solo_finder.merge_candidates(accepted, gap=0)[0]
    assert (candidate.bait, candidate.n_hits, candidate.end) == ("first", 2, 1499)


def test_hits_separated_by_more_than_the_gap_stay_separate() -> None:
    accepted = [
        solo_finder.Hit("b", "chr1", 99.0, 400, 400, 1000, 1399),
        solo_finder.Hit("b", "chr1", 99.0, 400, 400, 5000, 5399),
    ]
    assert len(solo_finder.merge_candidates(accepted, gap=0)) == 2


def test_candidates_on_different_sequences_never_merge() -> None:
    accepted = [
        solo_finder.Hit("b", "chr1", 99.0, 400, 400, 1000, 1399),
        solo_finder.Hit("b", "chr2", 99.0, 400, 400, 1000, 1399),
    ]
    assert len(solo_finder.merge_candidates(accepted, gap=0)) == 2


def test_reversed_subject_coordinates_are_normalised() -> None:
    """Blastn reports sstart > send for a minus-strand hit."""
    accepted = [solo_finder.Hit("b", "chr1", 99.0, 400, 400, 1399, 1000)]
    candidate = solo_finder.merge_candidates(accepted, gap=0)[0]
    assert (candidate.start, candidate.end) == (1000, 1399)


# ---- classification, and its priority order ----


def test_a_candidate_inside_an_intact_element_is_not_a_solo() -> None:
    elements = {"chr1": IntervalIndex([(49000, 61000)])}
    orphans = {"chr1": IntervalIndex([])}
    candidate = solo_finder.Candidate("chr1", 50000, 50399, 99.0, 1, "bait")
    assert solo_finder.classify(candidate, elements, orphans, orphan_pad=10000) == (
        solo_finder.INTACT_FLANK
    )


def test_a_candidate_near_an_orphan_is_a_mono_ltr_not_a_solo() -> None:
    elements = {"chr1": IntervalIndex([])}
    orphans = {"chr1": IntervalIndex([(95000, 95500)])}
    candidate = solo_finder.Candidate("chr1", 90000, 90399, 99.0, 1, "bait")
    assert solo_finder.classify(candidate, elements, orphans, orphan_pad=10000) == (
        solo_finder.MONO_AT_ORPHAN
    )


def test_a_candidate_that_is_neither_is_a_solo() -> None:
    elements = {"chr1": IntervalIndex([])}
    orphans = {"chr1": IntervalIndex([(95000, 95500)])}
    candidate = solo_finder.Candidate("chr1", 10000, 10399, 99.0, 1, "bait")
    assert (
        solo_finder.classify(candidate, elements, orphans, orphan_pad=10000)
        == solo_finder.SOLO
    )


def test_element_overlap_wins_over_orphan_proximity() -> None:
    """Both conditions true at once must resolve to intact, never to monoLTR.

    An intact element with a nearby orphan is a catalogued element; calling it a
    monoLTR would double-count it against the solo/intact ratio's denominator.
    """
    elements = {"chr1": IntervalIndex([(49000, 61000)])}
    orphans = {"chr1": IntervalIndex([(50100, 50200)])}
    candidate = solo_finder.Candidate("chr1", 50000, 50399, 99.0, 1, "bait")
    assert solo_finder.classify(candidate, elements, orphans, orphan_pad=10000) == (
        solo_finder.INTACT_FLANK
    )


def test_a_sequence_with_no_elements_or_orphans_yields_solos() -> None:
    """An unplaced scaffold carries neither track and must not crash the lookup."""
    candidate = solo_finder.Candidate("scaffold_9", 100, 499, 99.0, 1, "bait")
    assert solo_finder.classify(candidate, {}, {}, orphan_pad=10000) == solo_finder.SOLO


# ---- the whole funnel, end to end on the fixture ----


def test_the_funnel_counts_every_stage(
    hits: Path, elements: Path, orphans: Path, thresholds: solo_finder.Thresholds
) -> None:
    result = solo_finder.run(hits, elements, orphans, thresholds)
    funnel = result.funnel
    assert funnel["raw_hits"] == 5
    # Rejected: one on identity, one on alignment length, one on coverage.
    assert funnel["accepted_hits"] == 2
    # The two accepted hits are 10 kb apart, so they stay separate.
    assert funnel["merged_candidates"] == 2
    assert funnel[solo_finder.INTACT_FLANK] == 2
    assert funnel[solo_finder.SOLO] == 0


def test_rejection_reasons_are_tallied_for_diagnosis(
    hits: Path, elements: Path, orphans: Path, thresholds: solo_finder.Thresholds
) -> None:
    """The funnel plot needs to show where candidates die, not just how many."""
    result = solo_finder.run(hits, elements, orphans, thresholds)
    assert result.funnel["rejected_identity"] == 1
    assert result.funnel["rejected_alignment_length"] == 1
    assert result.funnel["rejected_coverage"] == 1


def test_an_empty_hit_table_is_a_valid_result_not_an_error(
    tmp_path: Path, elements: Path, orphans: Path, thresholds: solo_finder.Thresholds
) -> None:
    """Zero solos is legitimate for a genome with no ERVs."""
    empty = tmp_path / "empty.tsv.gz"
    with gzip.open(empty, "wt") as handle:
        handle.write("")
    result = solo_finder.run(empty, elements, orphans, thresholds)
    assert result.funnel["raw_hits"] == 0
    assert result.candidates == []


# ---- the solo list, which is the seam the integrator reads ----


def test_the_solo_list_has_the_six_columns_the_integrator_expects(
    tmp_path: Path,
) -> None:
    candidates = [
        solo_finder.ClassifiedCandidate(
            solo_finder.Candidate(
                "chr1", 1000, 1399, 99.0, 2, "chr1|LTR_retrotransposon1|L"
            ),
            solo_finder.SOLO,
            orphan_distance=50000,
        )
    ]
    spans = {"LTR_retrotransposon1": ("chr1", 49000, 61000)}
    out = tmp_path / "solo_list.tsv"
    solo_finder.write_solo_list(candidates, spans, out)
    fields = out.read_text().rstrip("\n").split("\t")
    assert len(fields) == 6
    assert fields[0] == "chr1"
    assert fields[1] == "1000"
    assert fields[2] == "1399"
    assert fields[3] == "chr1:1000..1399"
    # The library id carries the seeding element's span, which is how the
    # integrator maps a solo back to a classified locus by coordinate.
    assert fields[4] == "chr1:49000..61000#LTR/LTR_retrotransposon1"


def test_only_solos_reach_the_solo_list(tmp_path: Path) -> None:
    candidates = [
        solo_finder.ClassifiedCandidate(
            solo_finder.Candidate(
                "chr1", 1000, 1399, 99.0, 1, "chr1|LTR_retrotransposon1|L"
            ),
            solo_finder.INTACT_FLANK,
            orphan_distance=None,
        ),
        solo_finder.ClassifiedCandidate(
            solo_finder.Candidate(
                "chr1", 8000, 8399, 99.0, 1, "chr1|LTR_retrotransposon1|L"
            ),
            solo_finder.SOLO,
            orphan_distance=None,
        ),
    ]
    spans = {"LTR_retrotransposon1": ("chr1", 49000, 61000)}
    out = tmp_path / "solo_list.tsv"
    written = solo_finder.write_solo_list(candidates, spans, out)
    assert written == 1
    assert out.read_text().count("\n") == 1


def test_element_spans_are_read_from_the_elements_gff3(elements: Path) -> None:
    spans = solo_finder.element_spans(elements)
    assert spans == {"LTR_retrotransposon1": ("chr1", 49000, 61000)}
