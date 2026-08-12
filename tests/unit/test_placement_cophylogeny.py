"""Unit tests for ``workflow/scripts/taxonomy/placement_cophylogeny.py``.

Builds a tree of *host genomes* from their ERV placement distributions and sets
it against the host phylogeny. Congruent shapes mean ERVs were largely inherited
vertically with their hosts; discordant shapes point at cross-species
transmission or lineage-specific expansion and loss.

Two mechanics need pinning:

* **Sample naming.** gappa labels each sample by its file basename, and every
  EPA-ng run writes ``epa_result.jplace``. Passed directly, all five tips come
  back named ``epa_result`` and the tree is useless. Inputs are therefore staged
  under genome-derived names first.
* **Topology comparison.** Congruence is measured on bipartitions (the splits a
  tree induces), which is what Robinson-Foulds counts. Branch lengths are
  deliberately ignored: the supplied host tree is often a cladogram with
  placeholder lengths, and comparing those to placement distances would be
  meaningless.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from placement_cophylogeny import (
    bipartitions,
    congruence,
    sample_name_for,
    squash_cmd,
    stage_jplace,
)

# ((A,M),D),(H,Mu) - the host topology of the five model genomes
HOST = "(((Antrozous_pallidus:1,Molossus_molossus:1):1,Desmodus_rotundus:2):1,(Homo_sapiens:1,Mus_musculus:1):2);"
# What the ERV placements actually produced: Mus nested inside the bats
ERV = "(Homo_sapiens:0.65,(Molossus_molossus:0.46,(Mus_musculus:0.42,(Antrozous_pallidus:0.22,Desmodus_rotundus:0.22):0.21):0.15):0.16);"


# ---------------------------------------------------------------------
# sample naming
# ---------------------------------------------------------------------
def test_sample_name_is_the_genome_not_the_tier_or_gene() -> None:
    """Tips must be host genomes so the tree is comparable to the host phylogeny."""
    assert sample_name_for("Antrozous_pallidus.ltr-flanked.POL.jplace") == (
        "Antrozous_pallidus"
    )


def test_sample_name_handles_a_genome_containing_dots() -> None:
    """Splitting naively on '.' would truncate an assembly-style genome name."""
    assert sample_name_for("GCF_000001.2.orphan.POL.jplace") == "GCF_000001.2"


def test_stage_jplace_renames_inputs_to_genome_names(tmp_path: Path) -> None:
    """Without this every tip comes back as 'epa_result' (gappa uses basenames)."""
    src = tmp_path / "src"
    src.mkdir()
    for g in ("Antrozous_pallidus", "Mus_musculus"):
        (src / f"{g}.ltr-flanked.POL.jplace").write_text("{}")
    staged = stage_jplace(sorted(src.glob("*.jplace")), tmp_path / "staged")

    assert sorted(p.name for p in staged.iterdir()) == [
        "Antrozous_pallidus.jplace",
        "Mus_musculus.jplace",
    ]


def test_stage_jplace_refuses_duplicate_sample_names(tmp_path: Path) -> None:
    """Two tiers of one genome would collide and silently drop a sample.

    Tiers are analysed separately for exactly this reason; mixing them would
    compare a genome against itself.
    """
    src = tmp_path / "src"
    src.mkdir()
    (src / "Mus_musculus.ltr-flanked.POL.jplace").write_text("{}")
    (src / "Mus_musculus.orphan.POL.jplace").write_text("{}")
    with pytest.raises(ValueError, match="Mus_musculus"):
        stage_jplace(sorted(src.glob("*.jplace")), tmp_path / "staged")


# ---------------------------------------------------------------------
# command construction
# ---------------------------------------------------------------------
def test_squash_cmd_points_at_the_staged_directory(tmp_path: Path) -> None:
    cmd = squash_cmd(tmp_path / "staged", tmp_path / "out")
    assert cmd[:3] == ["gappa", "analyze", "squash"]
    assert "--write-newick-tree" in cmd
    assert str(tmp_path / "staged") in cmd


# ---------------------------------------------------------------------
# topology comparison
# ---------------------------------------------------------------------
def test_bipartitions_returns_the_nontrivial_splits() -> None:
    """A split is one side of an internal branch; trivial single-tip splits are
    excluded because every tree shares them and they carry no signal."""
    splits = bipartitions(HOST)
    assert frozenset({"Antrozous_pallidus", "Molossus_molossus"}) in splits
    assert frozenset({"Homo_sapiens", "Mus_musculus"}) in splits
    assert not any(len(s) < 2 for s in splits)


def test_identical_topologies_are_fully_congruent() -> None:
    result = congruence(HOST, HOST)
    assert result["shared_splits"] == result["host_splits"]
    assert result["rf_distance"] == 0
    assert result["congruent"] is True


def test_the_real_result_is_discordant() -> None:
    """The measured ERV tree nests Mus inside the bats, which the host tree does
    not - so this must not be reported as congruent."""
    result = congruence(HOST, ERV)
    assert result["congruent"] is False
    assert result["rf_distance"] > 0
    assert result["shared_splits"] < result["host_splits"]


def test_congruence_reports_which_grouping_disagrees() -> None:
    """A bare distance is not actionable; the caller needs the offending split."""
    result = congruence(HOST, ERV)
    text = " ".join(sorted(result["erv_only_splits"]))
    assert "Antrozous_pallidus" in text
    assert "Desmodus_rotundus" in text


def test_congruence_requires_a_shared_tip_set() -> None:
    """Comparing trees over different taxa silently produces nonsense splits."""
    with pytest.raises(ValueError, match="tip"):
        congruence(HOST, "(A:1,B:1);")
