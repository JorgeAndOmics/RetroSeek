"""Choose the tips for the solo-LTR evidence tree.

The tree answers one question the detector cannot: are the three fates of a lone
LTR real biological classes, or an artefact of how we drew the boundaries? If solos
were simply mis-called flanking arms they would scatter among them on an LTR
phylogeny. If some LTR families survive only as solos, solos will form their own
clades.

Tips are every bait arm plus a seeded sample of the other two classes. All arms are
kept because they are the reference frame and they carry the tree's positive
control: an element's two arms were identical the day it inserted, so they must come
out as sister tips. The other classes are sampled because a tree over every solo
would be neither computable nor readable, and the sample size and seed are config
values so the figure is reproducible.

Tip names encode the class as a prefix before a double underscore, which is what
`tree_stats.py` reads back:

    FLANK__{seqname}_{element}_{L|R}
    SOLO__{seqname}_{start}
    MONO__{seqname}_{start}
"""

from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path
from typing import NamedTuple

FLANK = "FLANK"
SOLO = "SOLO"
MONO = "MONO"

# Fate values as solo_finder writes them, mapped to the tree's class prefixes.
FATE_PREFIX = {"solo": SOLO, "mono_ltr_at_orphan": MONO}


class Tip(NamedTuple):
    seqname: str
    start: int  # BED, 0-based
    end: int
    name: str


def bait_tips(bait_bed: Path) -> list[Tip]:
    """Every bait arm, renamed with the FLANK prefix.

    The bait name is `{seqname}|{element}|{arm}`; `|` is replaced because a Newick
    label cannot carry it.
    """
    tips = []
    with bait_bed.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            tips.append(
                Tip(
                    seqname=fields[0],
                    start=int(fields[1]),
                    end=int(fields[2]),
                    name=f"{FLANK}__{fields[3].replace('|', '_')}",
                )
            )
    return tips


def sampled_tips(
    candidates_csv: Path, fate: str, n: int, rng: random.Random
) -> list[Tip]:
    """A seeded sample of candidates with the given fate.

    Sampling without replacement, and asking for more than exist yields all of them
    rather than raising: a genome can easily have fewer than `n` monoLTRs.
    """
    prefix = FATE_PREFIX[fate]
    pool: list[Tip] = []
    with candidates_csv.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row.get("fate") != fate:
                continue
            start, end = int(row["start"]), int(row["end"])
            pool.append(
                Tip(
                    seqname=row["seqname"],
                    start=start - 1,  # the candidate table is GFF3-style, 1-based
                    end=end,
                    name=f"{prefix}__{row['seqname']}_{start}",
                )
            )
    if n >= len(pool):
        return pool
    return rng.sample(pool, n)


def write_bed(tips: list[Tip], path: Path) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for tip in tips:
            handle.write(f"{tip.seqname}\t{tip.start}\t{tip.end}\t{tip.name}\t.\t+\n")
    return len(tips)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bait-bed", type=Path, required=True)
    parser.add_argument("--candidates-csv", type=Path, required=True)
    parser.add_argument("--out-bed", type=Path, required=True)
    parser.add_argument("--n-solo-tips", type=int, required=True)
    parser.add_argument("--n-mono-tips", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    args = parser.parse_args(argv)

    rng = random.Random(args.seed)
    tips = bait_tips(args.bait_bed)
    solos = sampled_tips(args.candidates_csv, "solo", args.n_solo_tips, rng)
    monos = sampled_tips(
        args.candidates_csv, "mono_ltr_at_orphan", args.n_mono_tips, rng
    )
    written = write_bed(tips + solos + monos, args.out_bed)
    print(
        f"tree tips: {len(tips)} flanking arms, {len(solos)} solos, "
        f"{len(monos)} monoLTRs-at-orphans -> {written} total"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
