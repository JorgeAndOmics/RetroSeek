"""Tests for solo-LTR tree tip selection.

The tree's positive control only means anything for an element whose BOTH arms are
on the tree, so sampling has to be by element, never by arm. And all three classes
have to be capped, because the clustering statistic compares against class
abundance: on Mus musculus, taking every bait arm made the tree 96.5% flanking,
raised the permutation null to 0.94, and collapsed the enrichment to 1.03x.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest
import solo_tips

_BAIT = "".join(
    f"chr1\t{i * 1000}\t{i * 1000 + 400}\tchr1|LTR_retrotransposon{i}|{arm}\t.\t+\n"
    for i in range(1, 11)
    for arm in ("L", "R")
)

_CANDIDATES = (
    "seqname,start,end,length,fate,best_identity,n_hits,bait,parent,orphan_distance\n"
    + "".join(
        f"chr2,{i * 500},{i * 500 + 400},400,{fate},99.0,1,chr1|LTR_retrotransposon1|L,"
        f"LTR_retrotransposon1,50000\n"
        for fate in ("solo", "mono_ltr_at_orphan")
        for i in range(1, 21)
    )
)


@pytest.fixture
def bait(tmp_path: Path) -> Path:
    p = tmp_path / "bait.bed"
    p.write_text(_BAIT)
    return p


@pytest.fixture
def candidates(tmp_path: Path) -> Path:
    p = tmp_path / "candidates.csv"
    p.write_text(_CANDIDATES)
    return p


def test_both_arms_of_a_sampled_element_are_kept(bait: Path) -> None:
    """Sampling by arm would break the control; sampling by element preserves it."""
    tips = solo_tips.bait_tips(bait, n_elements=3, rng=random.Random(1))
    elements = [t.name.rsplit("_", 1)[0] for t in tips]
    assert len(tips) == 6, "three elements must contribute two arms each"
    for element in set(elements):
        assert elements.count(element) == 2


def test_asking_for_more_elements_than_exist_keeps_them_all(bait: Path) -> None:
    tips = solo_tips.bait_tips(bait, n_elements=999, rng=random.Random(1))
    assert len(tips) == 20


def test_zero_means_no_cap(bait: Path) -> None:
    """0 is the documented way to keep every element."""
    tips = solo_tips.bait_tips(bait, n_elements=0, rng=random.Random(1))
    assert len(tips) == 20


def test_sampling_is_reproducible_for_a_seed(bait: Path) -> None:
    first = solo_tips.bait_tips(bait, n_elements=4, rng=random.Random(7))
    second = solo_tips.bait_tips(bait, n_elements=4, rng=random.Random(7))
    assert [t.name for t in first] == [t.name for t in second]


def test_tips_carry_their_class_prefix(bait: Path, candidates: Path) -> None:
    """tree_stats reads the class back off the prefix, so it has to be there."""
    rng = random.Random(1)
    flank = solo_tips.bait_tips(bait, n_elements=2, rng=rng)
    solos = solo_tips.sampled_tips(candidates, "solo", 5, rng)
    monos = solo_tips.sampled_tips(candidates, "mono_ltr_at_orphan", 5, rng)
    assert all(t.name.startswith("FLANK__") for t in flank)
    assert all(t.name.startswith("SOLO__") for t in solos)
    assert all(t.name.startswith("MONO__") for t in monos)


def test_sampled_candidates_are_converted_to_bed_coordinates(candidates: Path) -> None:
    """The candidate table is GFF3-style 1-based inclusive; BED starts one lower.

    The tip name keeps the ORIGINAL 1-based start, because that is what identifies
    the locus in the candidate table, so the two can be checked against each other.
    """
    tips = solo_tips.sampled_tips(candidates, "solo", 1, random.Random(1))
    tip = tips[0]
    locus = tip.name.split("__")[1]  # SOLO__{seqname}_{start}__{seed}
    gff_start = int(locus.rsplit("_", 1)[1])
    assert tip.start == gff_start - 1


def test_asking_for_more_candidates_than_exist_yields_all(candidates: Path) -> None:
    """A genome can easily have fewer monoLTRs than the cap."""
    tips = solo_tips.sampled_tips(
        candidates, "mono_ltr_at_orphan", 999, random.Random(1)
    )
    assert len(tips) == 20


# ---- the seed control (2026-09-22) ----
# Every solo was caught by one bait arm at >= 95% identity: its seed. A correct tree
# must put a solo beside its seed, which is the control that catches a broken tree
# (it caught mixed strands). So each solo names its seed, and seeds are always on
# the tree.


def test_a_solo_tip_names_its_seed_element(candidates: Path) -> None:
    tips = solo_tips.sampled_tips(candidates, "solo", 3, random.Random(1))
    assert all(t.name.endswith("__chr1_LTR_retrotransposon1") for t in tips)
    assert all(t.name.startswith("SOLO__chr2_") for t in tips)


def test_the_seeds_of_sampled_solos_are_always_on_the_tree(bait: Path) -> None:
    tips = solo_tips.bait_tips(
        bait, n_elements=2, rng=random.Random(1), required={"LTR_retrotransposon7"}
    )
    elements = {t.name.rsplit("_", 1)[0] for t in tips}
    assert "FLANK__chr1_LTR_retrotransposon7" in elements
    assert len(elements) == 2, "seeds count towards the cap; the rest is random"


def test_required_seeds_are_kept_even_beyond_the_cap(bait: Path) -> None:
    """The control needs every seed; the cap only limits the random fill."""
    required = {f"LTR_retrotransposon{i}" for i in (1, 2, 3)}
    tips = solo_tips.bait_tips(
        bait, n_elements=2, rng=random.Random(1), required=required
    )
    assert {
        t.name.rsplit("_", 1)[0].split("__")[1].split("_", 1)[1] for t in tips
    } == required


def test_a_solo_tip_carries_its_seed_element_id(candidates: Path) -> None:
    """A solo tip carries its seed element id as a field.

    It is not parsed back out of the name: sequence names such as RefSeq's NW_
    contain underscores, so the name cannot be split reliably.
    """
    tips = solo_tips.sampled_tips(candidates, "solo", 3, random.Random(1))
    assert {t.seed for t in tips} == {"LTR_retrotransposon1"}
