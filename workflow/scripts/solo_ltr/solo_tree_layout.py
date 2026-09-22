"""Lay out the solo-LTR evidence tree as coordinates, so R can draw it.

This is the same bridge ADR-011 built for the host and taxon trees: Python parses
the Newick and computes where every branch and tip goes, and writes two plain CSVs;
R reads them and draws segments and points. No R tree package is needed, and none
is in the environment, by the ADR's choice.

The layout itself is not reimplemented here. `taxonomy/tree_layout.layout` already
produces rectangular coordinates (x = distance from the root, y = ladderised tip
order) and is reused unchanged. What this adds is one column: each tip's class,
read from the prefix `solo_tips.py` put on its name (FLANK, SOLO or MONO), so the
figure can colour tips by fate.

Outputs:
    {genome}.tree_tips.csv      tip, class, x, y
    {genome}.tree_segments.csv  x, y, xend, yend
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import Phylo
from tree_layout import layout
from tree_stats import tip_class


def write_layout(treefile: Path, tips_csv: Path, segments_csv: Path) -> int:
    """Write the tip and segment coordinates; return the tip count.

    Branch lengths are kept (`align_tips=False`), because on an LTR tree they are
    the signal: a solo sitting on a long branch is a diverged copy, which is exactly
    what a reader wants to see.
    """
    tree = Phylo.read(str(treefile), "newick")  # type: ignore[no-untyped-call,attr-defined]
    segments, tips = layout(tree, align_tips=False)

    tips_csv.parent.mkdir(parents=True, exist_ok=True)
    with tips_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["tip", "class", "x", "y"])
        for name, x, y in tips:
            writer.writerow([name, tip_class(name or ""), x, y])

    segments_csv.parent.mkdir(parents=True, exist_ok=True)
    with segments_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["x", "y", "xend", "yend"])
        writer.writerows(segments)
    return len(tips)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--treefile", type=Path, required=True)
    parser.add_argument("--out-tips-csv", type=Path, required=True)
    parser.add_argument("--out-segments-csv", type=Path, required=True)
    args = parser.parse_args(argv)
    n = write_layout(args.treefile, args.out_tips_csv, args.out_segments_csv)
    print(f"tree layout: {n} tips -> {args.out_tips_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
