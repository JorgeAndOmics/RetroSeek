"""Cut the solo-LTR evidence tree into LTR families, and derive two views of it.

A family here is a maximal clade whose members are all close to each other: its
DIAMETER (the longest tip-to-tip path inside the clade, in substitutions per site)
is at most `max_distance`. Walking down from the root, the first clade that fits is
a family, and nothing inside it is split further.

Why 0.2 by default. It is the transposable-element convention: the 80-80-80 rule
(Wicker et al. 2007) groups sequences sharing at least 80% identity into one family,
and 0.2 substitutions per site is that boundary. Maximum-likelihood distances run
higher than raw mismatch, so the cut is slightly stricter than 80% identity. On the
pre-2026-09-23 trees it was also where family counts stopped tracking the cut; that
has not been re-measured on the corrected trees.

Each family gets a kind:

    no_intact    has solos and no sampled flanking arm. Since every sampled
                 solo's seed element is on the tree, this now arises only where
                 the tree has separated a solo from its seed: a diagnostic of
                 tree error, not a family without intact copies (every solo is
                 >= 95% identical to an intact element; see docs/solo_ltr.md).
    with_intact  has solos and at least one flanking arm
    no_solo      no solos at all

Two views are derived from the same tree, with no new inference:

    the solo-only tree   the tree pruned to its SOLO tips (the induced subtree)
    family subtrees      the largest `no_intact` and `with_intact` families, each
                         as its own small tree, so the two kinds can be compared

Tree depth is not a concern here: the model-5 trees are 38 to 47 nodes deep
(measured 2026-09-22), far inside Python's recursion limit that Biopython's own
traversals rely on.
"""

from __future__ import annotations

import copy
from collections import Counter
from dataclasses import dataclass
from typing import Any

from Bio.Phylo.BaseTree import Tree
from tree_stats import tip_class

NO_INTACT = "no_intact"
WITH_INTACT = "with_intact"
NO_SOLO = "no_solo"


@dataclass
class Family:
    """One family: its clade, identifier, and class census."""

    family_id: str
    clade: Any
    diameter: float
    n_flank: int
    n_solo: int
    n_mono: int

    @property
    def n_tips(self) -> int:
        """Tips of all three classes in the family."""
        return self.n_flank + self.n_solo + self.n_mono

    @property
    def kind(self) -> str:
        """no_solo without a solo; else with_intact or no_intact by flanking arms."""
        if self.n_solo == 0:
            return NO_SOLO
        return WITH_INTACT if self.n_flank else NO_INTACT


def annotate_diameters(tree: Any) -> None:
    """Give every clade a `height` (longest path down to a tip) and a `diameter`.

    Post-order, so children are always done before their parent. A clade's
    diameter is the larger of its children's diameters and the longest path that
    passes through the clade itself, which joins its two deepest children.
    """
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            clade.height, clade.diameter = 0.0, 0.0
            continue
        reach = sorted(
            (child.height + (child.branch_length or 0.0) for child in clade.clades),
            reverse=True,
        )
        clade.height = reach[0]
        through = reach[0] + reach[1] if len(reach) > 1 else reach[0]
        clade.diameter = max([through] + [child.diameter for child in clade.clades])


def cut_families(tree: Any, max_distance: float) -> list[Family]:
    """Split the tree into maximal clades of diameter <= max_distance.

    Families are numbered by size, largest first (ties broken by first tip name),
    so F001 is always the biggest and the numbering is stable across runs.
    """
    annotate_diameters(tree)
    clades = []
    stack = [tree.root]
    while stack:
        clade = stack.pop()
        if clade.diameter <= max_distance:
            clades.append(clade)
        else:
            stack.extend(clade.clades)

    def census(clade: Any) -> Counter[str]:
        return Counter(tip_class(tip.name or "") for tip in clade.get_terminals())

    counted = [(clade, census(clade)) for clade in clades]
    counted.sort(
        key=lambda item: (
            -sum(item[1].values()),
            min(t.name for t in item[0].get_terminals()),
        )
    )
    return [
        Family(
            family_id=f"F{i:03d}",
            clade=clade,
            diameter=clade.diameter,
            n_flank=counts["FLANK"],
            n_solo=counts["SOLO"],
            n_mono=counts["MONO"],
        )
        for i, (clade, counts) in enumerate(counted, 1)
    ]


def tip_families(families: list[Family]) -> dict[str, str]:
    """Tip name to family id."""
    return {tip.name: f.family_id for f in families for tip in f.clade.get_terminals()}


def solo_only_tree(tree: Any) -> Any | None:
    """The tree pruned to its SOLO tips; None when fewer than three remain.

    Pruning keeps the relationships the full tree inferred among these tips and
    merges branches where a pruned sibling leaves a single child, so path lengths
    between the remaining solos are unchanged.
    """
    pruned = copy.deepcopy(tree)
    for tip in list(pruned.get_terminals()):
        if tip_class(tip.name or "") != "SOLO":
            pruned.prune(tip)
    return pruned if pruned.count_terminals() >= 3 else None


def showcase(families: list[Family], per_kind: int) -> list[Family]:
    """The largest families of each solo-bearing kind, by solo count.

    Kept from the first design, which read no_intact families as families
    surviving only as solos; since the 2026-09-23 strand fix they flag tree error
    instead, and these pages await a redesign (backlog 9).
    """
    chosen = []
    for kind in (NO_INTACT, WITH_INTACT):
        ranked = sorted(
            (f for f in families if f.kind == kind and f.n_tips >= 3),
            key=lambda f: (-f.n_solo, f.family_id),
        )
        chosen.extend(ranked[:per_kind])
    return chosen


def as_tree(clade: Any) -> Any:
    """A standalone tree rooted at a copy of `clade`, for laying out on its own."""
    root = copy.deepcopy(clade)
    root.branch_length = 0.0
    return Tree(root=root, rooted=True)  # type: ignore[no-untyped-call]
