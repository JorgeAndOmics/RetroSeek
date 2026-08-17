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
    congruence_with_aliases,
    krd_cmd,
    output_prefix,
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
    cmd = squash_cmd(tmp_path / "staged", tmp_path / "out", "ltr-flanked.POL.")
    assert cmd[:3] == ["gappa", "analyze", "squash"]
    assert "--write-newick-tree" in cmd
    assert str(tmp_path / "staged") in cmd


# ---------------------------------------------------------------------
# Output naming
#
# gappa names its output files after the subcommand, not after the data:
# `analyze squash` always writes `cluster.newick` and `analyze krd` always
# writes `krd_matrix.csv`. Both tiers share one --out-dir, so with default
# names the two jobs write the same paths.
#
# Measured 2026-08-13: run concurrently, the tiers emitted a byte-identical
# composition tree (same md5) and both summaries reported the orphan splits.
# Run 27 minutes apart on 2026-08-12 they were correctly different. Snakemake
# could not catch it because the declared output, cophylogeny_summary.{tier}.
# {gene}.csv, *is* tier-scoped - only the undeclared intermediates collided.
# ---------------------------------------------------------------------
def test_output_prefix_distinguishes_tier_and_gene() -> None:
    """The prefix carries both, because both will vary independently.

    Only POL is placed today, but widening `placement_genes` is an open ADR-014
    follow-up, and a second gene would collide exactly as the tiers did.
    """
    assert output_prefix("ltr-flanked", "POL") == "ltr-flanked.POL."
    assert output_prefix("orphan", "POL") != output_prefix("ltr-flanked", "POL")
    assert output_prefix("orphan", "GAG") != output_prefix("orphan", "POL")


def test_squash_cmd_scopes_gappa_output_names_by_tier(tmp_path: Path) -> None:
    """The regression guard: two tiers must not name the same output file."""
    staged, out = tmp_path / "staged", tmp_path / "out"
    flanked = squash_cmd(staged, out, output_prefix("ltr-flanked", "POL"))
    orphan = squash_cmd(staged, out, output_prefix("orphan", "POL"))

    assert "--file-prefix" in flanked
    assert flanked[flanked.index("--file-prefix") + 1] == "ltr-flanked.POL."
    assert flanked != orphan


def test_krd_cmd_is_scoped_too(tmp_path: Path) -> None:
    """krd_matrix.csv collided the same way, leaving one tier's matrix only."""
    staged, out = tmp_path / "staged", tmp_path / "out"
    flanked = krd_cmd(staged, out, output_prefix("ltr-flanked", "POL"))
    orphan = krd_cmd(staged, out, output_prefix("orphan", "POL"))

    assert flanked[flanked.index("--file-prefix") + 1] == "ltr-flanked.POL."
    assert flanked != orphan


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


# ---------------------------------------------------------------------
# Alias-aware congruence
#
# The ERV-composition tree's tips are genome STEMS, taken from the published
# jplace filenames. A user's host tree is frequently labelled with display
# names, or with assembly directory names that match neither. Without
# canonicalisation the tip sets differ and the comparison is skipped entirely -
# reported as "not compared", which is honest but useless.
# ---------------------------------------------------------------------
SPECIES_MAP = {
    "GCF_014176215.1_mMyoMyo1": "Myotis myotis",
    "Antrozous_pallidus": "Antrozous pallidus",
    "GCA_004026805.1_ASM402680v1": "Desmodus rotundus",
}

ERV_STEMS = (
    "((GCF_014176215.1_mMyoMyo1:0.2,Antrozous_pallidus:0.2):0.1,"
    "GCA_004026805.1_ASM402680v1:0.3);"
)
# Display names in Newick must be quoted or underscored - unquoted whitespace is
# a token separator, so "Myotis myotis:12" parses as two names. Real exports use
# the underscore convention; both spellings are exercised below.
HOST_DISPLAY = (
    "(('Myotis myotis':12,'Antrozous pallidus':12):8,'Desmodus rotundus':20);"
)
HOST_UNDERSCORE = "((Myotis_myotis:12,Antrozous_pallidus:12):8,Desmodus_rotundus:20);"


def test_congruence_matches_display_names_against_genome_stems() -> None:
    """The real-workflow case: host tree in display names, ERV tree in stems."""
    result = congruence_with_aliases(HOST_DISPLAY, ERV_STEMS, SPECIES_MAP)
    assert result["n_tips"] == 3  # all three names resolved across spellings


def test_congruence_with_aliases_still_detects_a_real_disagreement() -> None:
    """Canonicalising names must not paper over a genuine topology difference.

    Uses four taxa: three tips admit only one unrooted topology, so any
    3-tip comparison is vacuously congruent.
    """
    smap = {**SPECIES_MAP, "Mus_musculus": "Mus musculus"}
    host = (
        "(('Myotis myotis':12,'Antrozous pallidus':12):8,"
        "('Desmodus rotundus':15,'Mus musculus':15):5);"
    )
    swapped = (
        "((GCF_014176215.1_mMyoMyo1:0.2,GCA_004026805.1_ASM402680v1:0.2):0.1,"
        "(Antrozous_pallidus:0.3,Mus_musculus:0.3):0.1);"
    )
    result = congruence_with_aliases(host, swapped, smap)
    assert result["comparable"] is True
    assert result["congruent"] is False
    assert result["rf_distance"] > 0


def test_three_taxa_are_reported_as_not_comparable() -> None:
    """Three tips have no internal branch, so agreement is vacuous.

    Reporting congruent=True there would look like a positive result when the
    comparison never actually happened.
    """
    result = congruence_with_aliases(HOST_UNDERSCORE, ERV_STEMS, SPECIES_MAP)
    assert result["comparable"] is False
    assert result["congruent"] is False


def test_congruence_matches_underscored_display_names_against_stems() -> None:
    """The common export convention: underscores standing in for spaces."""
    result = congruence_with_aliases(HOST_UNDERSCORE, ERV_STEMS, SPECIES_MAP)
    assert result["n_tips"] == 3  # underscored display names resolved too


def test_congruence_with_aliases_falls_back_when_no_map_is_given() -> None:
    """With no species map it delegates to the plain comparison unchanged."""
    assert congruence_with_aliases(HOST, ERV, {}) == congruence(HOST, ERV)
