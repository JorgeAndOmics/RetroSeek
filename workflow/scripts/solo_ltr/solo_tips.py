"""Choose the tips for the solo-LTR evidence tree.

The tree answers one question the detector cannot: are the three fates of a lone
LTR real biological classes, or an artefact of how we drew the boundaries? If solos
were simply mis-called flanking arms they would scatter among them on an LTR
phylogeny; if the fates are real, they cluster by LTR family, some families being
far richer in solos than others.

Tips are a seeded sample of all three classes. Flanking arms are the reference
frame and carry the tree's positive control (an element's two arms were identical
the day it inserted, so they tend to come out as sister tips), which is why they are
sampled by ELEMENT and both arms of a chosen element are kept.

All three classes are capped because the clustering statistic is a comparison
against class abundance, and it saturates when one class dominates. Taking every
bait arm did exactly that: on Mus musculus the tree came out 96.5% flanking arms,
the permutation null rose to 0.94, and the enrichment collapsed to 1.03x, which
looks like a result but is an artefact of composition. Capping elements keeps the
three classes within the same order of magnitude of each other and keeps the
control intact.

Tip names encode the class as a prefix before a double underscore, which is what
`tree_stats.py` reads back:

    FLANK__{seqname}_{element}_{L|R}
    SOLO__{seqname}_{start}__{seed seqname}_{seed element}
    MONO__{seqname}_{start}

A solo also names its SEED: the element whose LTR arm caught it at >= 95% identity.
Seeds are always sampled onto the tree, because a solo must sit beside its seed on
a correct tree. That is the tree's sharpest control: the arm control above could
not see the mixed-strand alignment of 2026-09 (an element's two arms share a
strand), while half the solos sat far from their seeds.
"""

from __future__ import annotations

import argparse
import csv
import logging
import random
from pathlib import Path
from typing import NamedTuple

from log import PipelineError, job_logging, run_main
from tabular import tab_rows

logger = logging.getLogger(__name__)

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
    seed: str = ""  # the seeding element's ID, for solos; empty otherwise


def _bait_arms(bait_bed: Path) -> dict[str, list[Tip]]:
    """Bait arms grouped by element, in file order.

    The bait name is `{seqname}|{element}|{arm}`; `|` is replaced in the tip name
    because a Newick label cannot carry it.
    """
    by_element: dict[str, list[Tip]] = {}
    with bait_bed.open() as handle:
        for fields in tab_rows(handle, 4):
            element = fields[3].split("|")[1] if "|" in fields[3] else fields[3]
            by_element.setdefault(element, []).append(
                Tip(
                    seqname=fields[0],
                    start=int(fields[1]),
                    end=int(fields[2]),
                    name=f"{FLANK}__{fields[3].replace('|', '_')}",
                )
            )
    return by_element


def bait_tips(
    bait_bed: Path,
    n_elements: int,
    rng: random.Random,
    required: set[str] | None = None,
) -> list[Tip]:
    """A seeded sample of bait arms, grouped so both arms of an element travel together.

    Sampling by element rather than by arm is what preserves the tree's positive
    control: it only means anything for an element whose two arms are both on the
    tree.

    `required` elements (the seeds of the sampled solos) are always kept, even past
    `n_elements`; the cap only limits the random fill around them.
    """
    by_element = _bait_arms(bait_bed)
    kept = sorted(set(required or ()) & set(by_element))
    others = sorted(set(by_element) - set(kept))
    room = len(others) if n_elements <= 0 else max(0, n_elements - len(kept))
    if room < len(others):
        others = rng.sample(others, room)
    return [tip for element in sorted(kept + others) for tip in by_element[element]]


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
            name = f"{prefix}__{row['seqname']}_{start}"
            seed = ""
            if prefix == SOLO:
                # bait is `{seqname}|{element}|{arm}`; the seed key drops the arm so
                # it matches the element part of FLANK__{seqname}_{element}_{arm}.
                seed = row["parent"]
                name += "__" + row["bait"].rsplit("|", 1)[0].replace("|", "_")
            pool.append(
                Tip(
                    seqname=row["seqname"],
                    start=start - 1,  # the candidate table is GFF3-style, 1-based
                    end=end,
                    name=name,
                    seed=seed,
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


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bait-bed", type=Path, required=True)
    parser.add_argument("--candidates-csv", type=Path, required=True)
    parser.add_argument("--out-bed", type=Path, required=True)
    parser.add_argument("--n-element-tips", type=int, required=True)
    parser.add_argument("--n-solo-tips", type=int, required=True)
    parser.add_argument("--n-mono-tips", type=int, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--log", type=Path, help="job log (the Snakemake log: path)")
    args = parser.parse_args(argv)
    job_logging(args.log, "solo_tree")

    rng = random.Random(args.seed)
    # Solos first, so their seeds can be required on the tree.
    solos = sampled_tips(args.candidates_csv, "solo", args.n_solo_tips, rng)
    monos = sampled_tips(
        args.candidates_csv, "mono_ltr_at_orphan", args.n_mono_tips, rng
    )
    seeds = {tip.seed for tip in solos}
    tips = bait_tips(args.bait_bed, args.n_element_tips, rng, required=seeds)
    written = write_bed(tips + solos + monos, args.out_bed)
    logger.info(
        "tree tips: %d flanking arms from at most %d elements (%d of them seeds of "
        "sampled solos), %d solos, %d monoLTRs at orphans; %d in all",
        len(tips),
        args.n_element_tips,
        len(seeds),
        len(solos),
        len(monos),
        written,
    )
    # IQ-TREE cannot build a tree from fewer than three sequences, and it fails
    # with a message that does not point back here. Say what is actually wrong.
    if written < 3:
        raise PipelineError(
            f"only {written} tree tips, too few to build a phylogeny: this genome "
            "has almost no ERV-bearing elements, so the tree would show nothing",
            hint="set solo_ltr.tree.enable to false to skip the tree stage",
        )


if __name__ == "__main__":
    run_main(main)
