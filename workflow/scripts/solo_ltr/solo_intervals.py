"""Interval overlap and distance lookup for the native solo-LTR detector.

The detector classifies each candidate locus by asking two questions:

1. does it overlap an LTRharvest element (then it is one flank of an intact
   element, not a solo);
2. how far is it from the nearest orphan locus (close means a monoLTR beside
   surviving retroviral coding sequence, so the provirus is damaged rather than
   excised).

Both are "does this interval touch that set of intervals", asked millions of times
against tens of thousands of intervals per chromosome, so the lookup has to be
sub-linear and it has to be exact.

Why this is its own module rather than a helper inside the detector: the prototype
got it wrong in a way that does not announce itself. It bisected the sorted starts
and then walked back at most 300 preceding intervals, which is correct only while
no more than 300 intervals start before the query. Past that it answers "no
overlap" for intervals that plainly overlap, and here a false "no overlap"
promotes an intact element to a solo LTR. Giving the lookup its own module and its
own tests makes that class of error visible.

The fix is to merge the input intervals on construction. Merged intervals are
disjoint and sorted, which makes both answers exact from a single bisect: if any
merged interval reaches the query, the last one starting at or before the query's
right edge is the one that does.
"""

from __future__ import annotations

import bisect
from collections.abc import Iterable


def merge_intervals(
    intervals: Iterable[tuple[int, int]], gap: int = 0
) -> list[tuple[int, int]]:
    """Sort and merge intervals, joining any pair no more than `gap` apart.

    `gap = 0` joins only intervals that overlap or touch. Coordinates are treated
    as inclusive on both ends, matching GFF3 rather than BED, because that is what
    both callers read.
    """
    ordered = sorted(intervals)
    if not ordered:
        return []
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + gap:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


class IntervalIndex:
    """Exact overlap and nearest-distance queries over one sequence's intervals.

    Construction merges, so `merged` is disjoint and ascending. Both queries are
    O(log n).
    """

    def __init__(self, intervals: Iterable[tuple[int, int]], gap: int = 0) -> None:
        self.merged = merge_intervals(intervals, gap=gap)
        self._starts = [start for start, _ in self.merged]

    def __len__(self) -> int:
        return len(self.merged)

    def overlaps(self, start: int, end: int, pad: int = 0) -> bool:
        """True if [start, end] expanded by `pad` touches any interval.

        Because the merged intervals are disjoint and ascending, the only one that
        can reach back to `start` is the last one beginning at or before
        `end + pad`: every earlier interval finishes before that one begins.
        """
        if not self.merged:
            return False
        i = bisect.bisect_right(self._starts, end + pad) - 1
        if i < 0:
            return False
        return self.merged[i][1] >= start - pad

    def nearest_distance(self, start: int, end: int) -> int | None:
        """Gap in bp to the nearest interval, 0 when overlapping, None when empty.

        Recorded as evidence rather than used as a filter (the ADR-015 principle):
        in a repeat-dense region a genuine solo can sit near an unrelated orphan by
        chance, so the distance belongs in a column where it can be inspected.
        """
        if not self.merged:
            return None
        i = bisect.bisect_right(self._starts, end) - 1
        best: int | None = None
        # The interval at or before the query, and the one after it. With disjoint
        # ascending intervals these are the only two candidates for nearest.
        for j in (i, i + 1):
            if not 0 <= j < len(self.merged):
                continue
            other_start, other_end = self.merged[j]
            if other_start <= end and other_end >= start:
                return 0
            gap = other_start - end if other_start > end else start - other_end
            best = gap if best is None else min(best, gap)
        return best
