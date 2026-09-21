"""Tests for the solo-LTR detector's interval index.

This module exists because of one specific defect. The prototype answered "does
this candidate touch any element" by bisecting the sorted starts and then walking
back at most 300 preceding intervals:

    j = bisect.bisect_right(starts, end + pad)
    for k in range(max(0, j - 300), j):   # <-- the bug
        ...

That is correct only while no more than 300 intervals start before the query. Past
that it silently answers "no overlap" for intervals that do overlap, and a silent
wrong answer here promotes an intact element to a solo LTR. The first test below is
the case that breaks it.
"""

from __future__ import annotations

import solo_intervals

# ---- the regression the module was written for ----


def test_overlap_is_found_beyond_the_prototypes_300_interval_window() -> None:
    """An overlap 400 entries back in the sorted list must still be found.

    The construction matters, so it is spelled out. One long interval, then 400
    short ones that all *start* before the query but *end* before it. They fill the
    backward-walk window without overlapping anything, so a 300-step walk stops at
    index 101 and never reaches index 0, which is the interval that does overlap.

    Verified to discriminate: the prototype's logic answers False here, and a
    candidate wrongly told it overlaps no element gets promoted to a solo LTR.
    """
    intervals = [(1, 10_000_000)]
    intervals += [
        (i * 100 + 10, i * 100 + 60) for i in range(400)
    ]  # all end below 40,000
    index = solo_intervals.IntervalIndex(intervals)

    assert index.overlaps(50_000, 51_000) is True


def test_no_overlap_is_still_no_overlap_with_many_intervals() -> None:
    """The fix must not turn every query into a hit."""
    intervals = [(i * 1000, i * 1000 + 100) for i in range(500)]
    index = solo_intervals.IntervalIndex(intervals)
    assert index.overlaps(500, 900) is False


# ---- boundaries ----


def test_touching_at_a_single_base_counts_as_overlap() -> None:
    index = solo_intervals.IntervalIndex([(100, 200)])
    assert index.overlaps(200, 300) is True
    assert index.overlaps(50, 100) is True


def test_adjacent_but_not_touching_does_not_overlap() -> None:
    index = solo_intervals.IntervalIndex([(100, 200)])
    assert index.overlaps(201, 300) is False
    assert index.overlaps(50, 99) is False


def test_a_nested_interval_is_found() -> None:
    """Nesting is why the index merges on construction rather than trusting order."""
    index = solo_intervals.IntervalIndex([(100, 10_000), (200, 300)])
    assert index.overlaps(5_000, 5_001) is True


def test_pad_extends_the_reach_symmetrically() -> None:
    index = solo_intervals.IntervalIndex([(1_000, 2_000)])
    assert index.overlaps(2_500, 2_600) is False
    assert index.overlaps(2_500, 2_600, pad=500) is True
    assert index.overlaps(400, 500, pad=500) is True


def test_an_empty_index_never_overlaps() -> None:
    index = solo_intervals.IntervalIndex([])
    assert index.overlaps(1, 100) is False
    assert index.nearest_distance(1, 100) is None


# ---- distance, which is recorded as evidence rather than used as a gate ----


def test_distance_is_zero_when_overlapping() -> None:
    index = solo_intervals.IntervalIndex([(100, 200)])
    assert index.nearest_distance(150, 160) == 0


def test_distance_looks_both_ways_and_takes_the_smaller_gap() -> None:
    index = solo_intervals.IntervalIndex([(100, 200), (1_000, 1_100)])
    # 300..400 is 100 from the left interval and 600 from the right one.
    assert index.nearest_distance(300, 400) == 100
    # 900..950 is 700 from the left and 50 from the right.
    assert index.nearest_distance(900, 950) == 50


def test_distance_past_the_last_interval_is_measured_from_its_end() -> None:
    index = solo_intervals.IntervalIndex([(100, 200)])
    assert index.nearest_distance(1_000, 1_100) == 800


# ---- merging ----


def test_overlapping_inputs_are_merged() -> None:
    index = solo_intervals.IntervalIndex([(100, 200), (150, 300), (290, 400)])
    assert index.merged == [(100, 400)]


def test_touching_inputs_are_merged() -> None:
    index = solo_intervals.IntervalIndex([(100, 200), (200, 300)])
    assert index.merged == [(100, 300)]


def test_disjoint_inputs_are_kept_apart_and_sorted() -> None:
    index = solo_intervals.IntervalIndex([(500, 600), (100, 200)])
    assert index.merged == [(100, 200), (500, 600)]
