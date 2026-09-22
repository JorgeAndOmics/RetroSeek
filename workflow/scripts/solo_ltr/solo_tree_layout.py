"""Lay out the solo-LTR evidence tree as coordinates, so R can draw it.

This is the same bridge ADR-011 built for the host and taxon trees: Python parses
the Newick and computes where every branch and tip goes, and writes plain CSVs; R
reads them and draws segments and points. No R tree package is needed, and none is
in the environment, by the ADR's choice. The rectangular layout itself is
`taxonomy/tree_layout.layout`, reused unchanged.

Three views of the one tree are written, none needing a new inference:

    the full tree        every tip, coloured by fate in the figure
    the solo-only tree   the tree pruned to its solos, coloured by family
    family subtrees      the largest no-intact and with-intact families, each
                         laid out on its own (see tree_families.py)

plus the family table that ties them together.

Outputs, for `--table-prefix {dir}/{genome}`:
    {genome}.tree_tips.csv              tip, class, family, x, y
    {genome}.tree_segments.csv          x, y, xend, yend
    {genome}.solo_tree_tips.csv         tip, family, x, y
    {genome}.solo_tree_segments.csv     x, y, xend, yend
    {genome}.family_tree_tips.csv       family, tip, class, x, y
    {genome}.family_tree_segments.csv   family, x, y, xend, yend
    {genome}.tree_families.csv          one row per family, with its kind
and the pruned tree as Newick at `--out-solo-newick`, for a tree viewer.
"""

from __future__ import annotations

import argparse
import csv
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import tree_families
from Bio import Phylo
from tree_layout import layout
from tree_stats import tip_class


def _write(path: Path, header: list[str], rows: Iterable[Iterable[Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def write_views(
    treefile: Path,
    table_prefix: Path,
    solo_newick: Path,
    max_distance: float,
    per_kind: int,
) -> list[tree_families.Family]:
    """Write every view and return the families, for the caller's summary line.

    Branch lengths are kept (`align_tips=False`) in every view, because on an LTR
    tree they are the signal: a solo on a long branch is a diverged copy.
    """
    tree = Phylo.read(str(treefile), "newick")  # type: ignore[no-untyped-call,attr-defined]
    families = tree_families.cut_families(tree, max_distance)
    family_of = tree_families.tip_families(families)

    def out(suffix: str) -> Path:
        return table_prefix.with_name(f"{table_prefix.name}.{suffix}")

    # The full tree.
    segments, tips = layout(tree, align_tips=False)
    _write(
        out("tree_tips.csv"),
        ["tip", "class", "family", "x", "y"],
        ((name, tip_class(name), family_of[name], x, y) for name, x, y in tips),
    )
    _write(out("tree_segments.csv"), ["x", "y", "xend", "yend"], segments)

    # The solo-only tree. Header-only files when there are too few solos, so the
    # outputs always exist and the figure shows a placeholder instead of failing.
    solo_tree = tree_families.solo_only_tree(tree)
    solo_segments: list[Any] = []
    solo_tips: list[Any] = []
    if solo_tree is not None:
        solo_segments, solo_tips = layout(solo_tree, align_tips=False)
        Phylo.write(solo_tree, str(solo_newick), "newick")  # type: ignore[no-untyped-call,attr-defined]
    else:
        solo_newick.write_text("")
    _write(
        out("solo_tree_tips.csv"),
        ["tip", "family", "x", "y"],
        ((name, family_of[name], x, y) for name, x, y in solo_tips),
    )
    _write(out("solo_tree_segments.csv"), ["x", "y", "xend", "yend"], solo_segments)

    # The showcase families, each laid out as its own tree.
    shown = tree_families.showcase(families, per_kind)
    family_tips: list[Any] = []
    family_segments: list[Any] = []
    for family in shown:
        segs, ftips = layout(tree_families.as_tree(family.clade), align_tips=False)
        family_tips.extend(
            (family.family_id, n, tip_class(n), x, y) for n, x, y in ftips
        )
        family_segments.extend((family.family_id, *seg) for seg in segs)
    _write(
        out("family_tree_tips.csv"), ["family", "tip", "class", "x", "y"], family_tips
    )
    _write(
        out("family_tree_segments.csv"),
        ["family", "x", "y", "xend", "yend"],
        family_segments,
    )

    shown_ids = {f.family_id for f in shown}
    _write(
        out("tree_families.csv"),
        [
            "family",
            "kind",
            "n_tips",
            "n_flank",
            "n_solo",
            "n_mono",
            "diameter",
            "shown",
        ],
        (
            (
                f.family_id,
                f.kind,
                f.n_tips,
                f.n_flank,
                f.n_solo,
                f.n_mono,
                round(f.diameter, 4),
                f.family_id in shown_ids,
            )
            for f in families
        ),
    )
    return families


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--treefile", type=Path, required=True)
    parser.add_argument("--table-prefix", type=Path, required=True)
    parser.add_argument("--out-solo-newick", type=Path, required=True)
    parser.add_argument("--family-max-distance", type=float, required=True)
    parser.add_argument("--family-panels-per-kind", type=int, required=True)
    args = parser.parse_args(argv)

    families = write_views(
        args.treefile,
        args.table_prefix,
        args.out_solo_newick,
        args.family_max_distance,
        args.family_panels_per_kind,
    )
    no_intact = [f for f in families if f.kind == tree_families.NO_INTACT]
    with_intact = [f for f in families if f.kind == tree_families.WITH_INTACT]
    total_solos = sum(f.n_solo for f in families)
    print(
        f"families at diameter <= {args.family_max_distance}: {len(families)}; "
        f"{len(no_intact)} without an intact member, holding {sum(f.n_solo for f in no_intact)} of "
        f"{total_solos} sampled solos, {len(with_intact)} with an intact member"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
