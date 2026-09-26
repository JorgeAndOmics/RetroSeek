"""Measure what the solo-LTR tree says, rather than eyeballing it.

Three numbers, each answering a question a reader would otherwise have to take on
trust.

**The arm control.** An element's two LTR arms were identical the day it inserted,
so they tend to come out as sister tips. The fraction recovered as sisters says
whether the tree carries signal, but it is blurred wherever a family burst produced
many near-identical copies: an arm's sister is then as likely another element's arm
(about 61% on the model 5). It is also blind to strand errors, since an element's
two arms always share a strand.

**Do the classes cluster?** For each tip, does its sister group contain at least one
tip of its own class? Compared against a null built by permuting the class labels
over the fixed topology, which is what turns "solos seem to group together" into a
measurement. Without the null the statistic is meaningless, because a tree with any
structure at all will show some clustering.

**The seed control.** Every solo was caught by one bait arm at >= 95% identity, its
seed, so on a correct tree a solo sits within about 0.05 substitutions/site of its
seed's arms. The fraction within `SEED_CONTROL_DISTANCE` is the sharpest check the
tree gets: it exposed the mixed-strand alignment of 2026-09 (half the solos sat 0.2
to 3.6 away), which the arm control could not see, since an element's two arms
always share a strand.

**Who sits next to whom.** The class-by-class adjacency, as enrichment over what
class abundance alone predicts. Solos enriched beside other solos says the fates are
not spread evenly over LTR families: some families hold many solos per intact
element. It cannot mean families surviving ONLY as solos: every solo matches its
seed, an intact ERV-bearing element, at >= 95% identity, so none is without an
intact relative in the genome.

The tree itself is exploratory: `-fast`, no bootstrap, and a sample of the solos. It
is evidence that the classes are real, never an input to classification.
"""

from __future__ import annotations

import argparse
import csv
import logging
import random
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from Bio import Phylo

from log import OK, job_logging, run_main

logger = logging.getLogger(__name__)

# FLANK tips are named FLANK__{seqname}_{element}_{L|R}; the element is whatever
# sits between the prefix and the final arm letter.
ARM_PATTERN = re.compile(r"^FLANK__(.+)_(L|R)$")
CLASSES = ("FLANK", "SOLO", "MONO")

# A solo within this many substitutions/site of its seed's nearer arm passes the
# seed control. 95% identity is about 0.05; twice that absorbs alignment and
# model noise, and is still far below the 0.2 family cut.
SEED_CONTROL_DISTANCE = 0.1
# Below this share of solos near their seed, the tree is not to be trusted: the
# model 5 sit at 0.92 to 1.00 once sequences are oriented (2026-09-23).
SEED_CONTROL_MIN_FRACTION = 0.9


def tip_class(name: str) -> str:
    """The class prefix a tip carries, e.g. FLANK__chr1_x_L -> FLANK."""
    return name.split("__", 1)[0]


def parent_map(tree: Any) -> dict[Any, Any]:
    """Child clade to parent clade. Bio.Phylo has no upward links of its own."""
    parents = {}
    for clade in tree.find_clades(order="level"):
        for child in clade.clades:
            parents[child] = clade
    return parents


def sister_leaves(clade: Any, parents: dict[Any, Any]) -> list[Any]:
    """The terminals of every sibling subtree of `clade`."""
    parent = parents.get(clade)
    if parent is None:
        return []
    leaves: list[Any] = []
    for sibling in parent.clades:
        if sibling is clade:
            continue
        leaves.extend([sibling] if sibling.is_terminal() else sibling.get_terminals())
    return leaves


def arm_sisterhood(tree: Any, parents: dict[Any, Any]) -> tuple[int, int]:
    """(elements with both arms on the tree, of those recovered as sister tips)."""
    by_element: dict[str, list[Any]] = defaultdict(list)
    for tip in tree.get_terminals():
        match = ARM_PATTERN.match(tip.name or "")
        if match:
            by_element[match.group(1)].append(tip)
    both = [arms for arms in by_element.values() if len(arms) == 2]
    sisters = sum(1 for a, b in both if parents.get(a) is parents.get(b))
    return len(both), sisters


def seed_distances(tree: Any) -> list[float]:
    """Tree distance from each solo to the nearer arm of its seed, where on the tree.

    Solo tips are named SOLO__{locus}__{seed}, the seed key matching the element
    part of the seed's FLANK__{seed}_{L|R} arm names. Solos without a seed on the
    tree (or from before seeds were named) are skipped.
    """
    arms: dict[str, list[Any]] = defaultdict(list)
    for tip in tree.get_terminals():
        match = ARM_PATTERN.match(tip.name or "")
        if match:
            arms[match.group(1)].append(tip)
    distances = []
    for tip in tree.get_terminals():
        parts = (tip.name or "").split("__")
        if parts[0] != "SOLO" or len(parts) < 3 or parts[2] not in arms:
            continue
        distances.append(min(tree.distance(tip, arm) for arm in arms[parts[2]]))
    return distances


def same_class_sister_fraction(
    tips: list[Any], parents: dict[Any, Any], labels: dict[str, str]
) -> float:
    """Fraction of tips whose sister group contains their own class."""
    hits = considered = 0
    for tip in tips:
        leaves = sister_leaves(tip, parents)
        if not leaves:
            continue
        considered += 1
        if any(labels[leaf.name] == labels[tip.name] for leaf in leaves):
            hits += 1
    return hits / considered if considered else 0.0


def permutation_null(
    tips: list[Any],
    parents: dict[Any, Any],
    labels: dict[str, str],
    permutations: int,
    seed: int,
) -> list[float]:
    """The same statistic with class labels shuffled over the fixed topology."""
    rng = random.Random(seed)
    names = list(labels)
    values = [labels[name] for name in names]
    null = []
    for _ in range(permutations):
        rng.shuffle(values)
        null.append(
            same_class_sister_fraction(
                tips, parents, dict(zip(names, values, strict=True))
            )
        )
    return null


def _sister_pairs(tips: list[Any], parents: dict[Any, Any]) -> dict[str, Counter[str]]:
    """For each tip class, how often each class sits among its sisters."""
    pairs: dict[str, Counter[str]] = defaultdict(Counter)
    for tip in tips:
        own = tip_class(tip.name)
        for leaf in sister_leaves(tip, parents):
            pairs[own][tip_class(leaf.name)] += 1
    return pairs


def _adjacency_row(
    own: str, other: str, observed: Counter[str], total: int, expected_fraction: float
) -> dict[str, Any]:
    """One (tip class, sister class) row; enrichment is blank when not expected.

    ``total`` is the number of sisters seen for ``own``, over every class.
    """
    observed_fraction = (observed[other] / total) if total else 0.0
    return {
        "tip_class": own,
        "sister_class": other,
        "n": observed[other],
        "observed_fraction": round(observed_fraction, 4),
        "expected_fraction": round(expected_fraction, 4),
        "enrichment": (
            round(observed_fraction / expected_fraction, 3) if expected_fraction else ""
        ),
    }


def adjacency(tips: list[Any], parents: dict[Any, Any]) -> list[dict[str, Any]]:
    """Class-by-class sister counts, with enrichment over class abundance.

    The expectation is each class's share of all tips: if adjacency were random, a
    tip's sisters would be drawn in proportion to how common each class is.
    """
    pairs = _sister_pairs(tips, parents)
    abundance = Counter(tip_class(tip.name) for tip in tips)
    total_tips = sum(abundance.values()) or 1
    totals = {own: sum(pairs[own].values()) for own in CLASSES}
    return [
        _adjacency_row(
            own, other, pairs[own], totals[own], abundance[other] / total_tips
        )
        for own in CLASSES
        for other in CLASSES
    ]


def per_class_sisterhood(
    tips: list[Any], parents: dict[Any, Any]
) -> list[dict[str, Any]]:
    """Same-class-sister rate broken down by class."""
    counts: dict[str, list[int]] = defaultdict(lambda: [0, 0])
    for tip in tips:
        sisters = [tip_class(leaf.name) for leaf in sister_leaves(tip, parents)]
        if not sisters:
            continue
        own = tip_class(tip.name)
        counts[own][1] += 1
        if own in sisters:
            counts[own][0] += 1
    return [
        {
            "tip_class": name,
            "tips": counts[name][1],
            "same_class_sister": counts[name][0],
            "fraction": round(counts[name][0] / counts[name][1], 4)
            if counts[name][1]
            else "",
        }
        for name in CLASSES
        if name in counts
    ]


def summarise(treefile: Path, permutations: int, seed: int) -> dict[str, Any]:
    """Every statistic, plus the tip census, from one Newick file."""
    tree = Phylo.read(str(treefile), "newick")  # type: ignore[no-untyped-call,attr-defined]
    tips = tree.get_terminals()
    parents = parent_map(tree)
    labels = {tip.name: tip_class(tip.name) for tip in tips}

    both_arms, sisters = arm_sisterhood(tree, parents)
    observed = same_class_sister_fraction(tips, parents, labels)
    null = permutation_null(tips, parents, labels, permutations, seed)
    null_mean = statistics.mean(null) if null else 0.0
    null_sd = statistics.pstdev(null) if len(null) > 1 else 0.0

    census = Counter(labels.values())
    to_seed = seed_distances(tree)
    return {
        "n_tips": len(tips),
        "n_flank": census["FLANK"],
        "n_solo": census["SOLO"],
        "n_mono": census["MONO"],
        "elements_with_both_arms": both_arms,
        "arms_recovered_as_sisters": sisters,
        "arm_sisterhood_fraction": round(sisters / both_arms, 4) if both_arms else "",
        "solos_with_seed_on_tree": len(to_seed),
        "solos_near_seed_fraction": (
            round(sum(d <= SEED_CONTROL_DISTANCE for d in to_seed) / len(to_seed), 4)
            if to_seed
            else ""
        ),
        "seed_distance_median": round(statistics.median(to_seed), 4) if to_seed else "",
        "same_class_sister_observed": round(observed, 4),
        "same_class_sister_null_mean": round(null_mean, 4),
        "same_class_sister_null_sd": round(null_sd, 4),
        "enrichment": round(observed / null_mean, 3) if null_mean else "",
        "sd_above_null": round((observed - null_mean) / null_sd, 1) if null_sd else "",
        "permutations": permutations,
        "seed": seed,
        "adjacency": adjacency(tips, parents),
        "per_class": per_class_sisterhood(tips, parents),
    }


def write_csvs(summary: dict[str, Any], summary_csv: Path, adjacency_csv: Path) -> None:
    """Long-format summary plus the adjacency table, as the plots read them."""
    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    with summary_csv.open("w", newline="") as handle:
        summary_writer = csv.writer(handle)
        summary_writer.writerow(["metric", "value"])
        for key, value in summary.items():
            if key in {"adjacency", "per_class"}:
                continue
            summary_writer.writerow([key, value])
        for row in summary["per_class"]:
            summary_writer.writerow(
                [f"same_class_sister_{row['tip_class']}", row["fraction"]]
            )

    adjacency_csv.parent.mkdir(parents=True, exist_ok=True)
    with adjacency_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "tip_class",
                "sister_class",
                "n",
                "observed_fraction",
                "expected_fraction",
                "enrichment",
            ],
        )
        writer.writeheader()
        writer.writerows(summary["adjacency"])


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--treefile", type=Path, required=True)
    parser.add_argument("--out-summary-csv", type=Path, required=True)
    parser.add_argument("--out-adjacency-csv", type=Path, required=True)
    parser.add_argument("--permutations", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--log", type=Path, help="job log (the Snakemake log: path)")
    args = parser.parse_args(argv)
    job_logging(args.log, "solo_tree")

    summary = summarise(args.treefile, args.permutations, args.seed)
    write_csvs(summary, args.out_summary_csv, args.out_adjacency_csv)
    logger.info(
        "tree: %s tips (FLANK %s, SOLO %s, MONO %s)",
        summary["n_tips"],
        summary["n_flank"],
        summary["n_solo"],
        summary["n_mono"],
    )
    logger.info(
        "controls: arms recovered as sisters %s of %s (%s); %s of %s solos within "
        "%s of their seed",
        summary["arms_recovered_as_sisters"],
        summary["elements_with_both_arms"],
        summary["arm_sisterhood_fraction"],
        summary["solos_near_seed_fraction"],
        summary["solos_with_seed_on_tree"],
        SEED_CONTROL_DISTANCE,
    )
    logger.info(
        "clustering: observed %s against a null of %s (sd %s), %sx, %s sd above it",
        summary["same_class_sister_observed"],
        summary["same_class_sister_null_mean"],
        summary["same_class_sister_null_sd"],
        summary["enrichment"],
        summary["sd_above_null"],
    )
    near_seed = summary["solos_near_seed_fraction"]
    if isinstance(near_seed, float) and near_seed < SEED_CONTROL_MIN_FRACTION:
        logger.warning(
            "only %.0f%% of solos sit within %s of their seed arm (expected at least "
            "%.0f%%), so this genome's evidence tree is unreliable. Read its tree "
            "pages with care; a slower tree search (solo_ltr.tree.fast: false) may "
            "help",
            100 * near_seed,
            SEED_CONTROL_DISTANCE,
            100 * SEED_CONTROL_MIN_FRACTION,
        )
    logger.log(
        OK,
        "evidence tree: %s tips, %s of solos near their seed",
        summary["n_tips"],
        near_seed,
    )


if __name__ == "__main__":
    run_main(main)
