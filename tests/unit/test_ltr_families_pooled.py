"""Tests for pooled LTR families across genomes (solo_ltr/ltr_families_pooled.py).

The same clustering as the per-genome families, over every genome's bait arms at
once, so a family shared by two bats is one family. What can go wrong quietly: an
arm name reused by two genomes (chromosome names can repeat), and identifiers that
depend on the order the genomes are listed in.
"""

from __future__ import annotations

from pathlib import Path

import ltr_families_pooled as pooled
import pytest
from ltr_families import parse_clstr


def _bait(
    folder: Path, genome: str, arms: dict[str, tuple[str, int, int, str]]
) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f"{genome}.bait.fna").write_text(
        "".join(f">{name}\n{seq}\n" for name, (_, _, _, seq) in arms.items())
    )
    (folder / f"{genome}.bait.bed").write_text(
        "".join(
            f"{seqname}\t{start - 1}\t{end}\t{name}\t.\t+\n"
            for name, (seqname, start, end, _) in arms.items()
        )
    )


@pytest.fixture
def bait_dir(tmp_path: Path) -> Path:
    folder = tmp_path / "bait"
    _bait(
        folder,
        "Mus_musculus",
        {"chr1|LTR_retrotransposon1|L": ("chr1", 100, 400, "ACGT")},
    )
    _bait(
        folder,
        "Homo_sapiens",
        {
            "chr1|LTR_retrotransposon1|L": ("chr1", 50, 350, "GGCC"),
            "chr2|LTR_retrotransposon7|R": ("chr2", 10, 300, "TTAA"),
        },
    )
    return folder


def test_arms_of_every_genome_are_pooled_under_prefixed_names(
    bait_dir: Path, tmp_path: Path
) -> None:
    out = tmp_path / "pooled.fna"
    arms = pooled.pool_bait(bait_dir, ["Mus_musculus", "Homo_sapiens"], out)
    headers = [
        line[1:] for line in out.read_text().splitlines() if line.startswith(">")
    ]
    # Genomes in sorted order: cd-hit's result can depend on input order.
    assert headers == [
        "Homo_sapiens|chr1|LTR_retrotransposon1|L",
        "Homo_sapiens|chr2|LTR_retrotransposon7|R",
        "Mus_musculus|chr1|LTR_retrotransposon1|L",
    ]
    # The same arm name in two genomes stays two arms.
    assert len(arms) == 3
    assert arms["Homo_sapiens|chr1|LTR_retrotransposon1|L"].start == 50


def test_pooled_families_name_their_genomes(bait_dir: Path, tmp_path: Path) -> None:
    arms = pooled.pool_bait(
        bait_dir, ["Mus_musculus", "Homo_sapiens"], tmp_path / "p.fna"
    )
    clstr = tmp_path / "p.clstr"
    clstr.write_text(
        ">Cluster 0\n"
        "0\t300nt, >Homo_sapiens|chr1|LTR_retrotransposon1|L... *\n"
        "1\t300nt, >Mus_musculus|chr1|LTR_retrotransposon1|L... at +/85.00%\n"
        ">Cluster 1\n"
        "0\t290nt, >Homo_sapiens|chr2|LTR_retrotransposon7|R... *\n"
    )
    families = pooled.name_pooled(parse_clstr(clstr), arms)
    rows = pooled.pooled_rows(families)
    assert rows[0] == {
        "genome": "Homo_sapiens",
        "arm": "chr1|LTR_retrotransposon1|L",
        "pool_family": "Pool_F001",
        "representative": True,
    }
    summary = pooled.pooled_summary(families)
    assert summary[0]["n_genomes"] == 2
    assert summary[0]["arms_per_genome"] == "Homo_sapiens:1;Mus_musculus:1"
    assert summary[1]["pool_family"] == "Pool_F002"


def test_identifiers_do_not_depend_on_genome_order(
    bait_dir: Path, tmp_path: Path
) -> None:
    one = pooled.pool_bait(
        bait_dir, ["Mus_musculus", "Homo_sapiens"], tmp_path / "a.fna"
    )
    two = pooled.pool_bait(
        bait_dir, ["Homo_sapiens", "Mus_musculus"], tmp_path / "b.fna"
    )
    assert one == two
    assert (tmp_path / "a.fna").read_text() == (tmp_path / "b.fna").read_text()


def test_two_genomes_with_one_family_code_still_pool(tmp_path: Path) -> None:
    # Canis_lupus_familiaris and Canis_lupus_dingo both give Clup; the pooled step
    # names arms by the full genome name, so the codes never meet here.
    folder = tmp_path / "bait"
    arm = {"chr1|LTR_retrotransposon1|L": ("chr1", 1, 4, "ACGT")}
    _bait(folder, "Canis_lupus_familiaris", arm)
    _bait(folder, "Canis_lupus_dingo", arm)
    arms = pooled.pool_bait(
        folder, ["Canis_lupus_familiaris", "Canis_lupus_dingo"], tmp_path / "p.fna"
    )
    assert len(arms) == 2
