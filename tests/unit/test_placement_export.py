"""Unit tests for placement-artifact export (``taxonomy_placement.export_placement``).

EPA-ng writes ``epa_result.jplace`` and gappa writes ``labelled_tree.newick``
into a scratch workdir under ``data/tmp/``, documented as "cleared between runs".
That jplace is the evidence behind every ``taxon_call`` in the catalog - which
branch of the retroviral phylogeny each locus attached to, with likelihood
weights - and it is also the standard interchange format iTOL and gappa consume.
Discarding it means a reviewer cannot be shown why a locus was called.

Export copies those artifacts somewhere durable. The awkward case is that
placement legitimately does not run in three situations: the gene has no tree
package, every query aligned to all-gaps, or no locus carried that gene at all.
A Snakemake rule that declares the jplace as an output still needs a file in all
three, so export synthesises an empty-but-valid one carrying the real reference
tree. `tree_layout.py` sets the same precedent by writing header-only CSVs so
the DAG holds when no tree is configured.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest
from taxonomy_placement import edge_numbered_newick, empty_jplace, export_placement

REF_TREE = "((A:0.1,B:0.2):0.3,(C:0.15,D:0.25):0.35);"


def _write_ref_package(ref_dir: Path, gene: str = "POL") -> Path:
    """Minimal reference tree package: the two files export looks for."""
    trees = ref_dir / "trees"
    trees.mkdir(parents=True, exist_ok=True)
    (trees / f"{gene}.treefile").write_text(REF_TREE + "\n")
    (trees / f"{gene}.raxml.bestTree").write_text(REF_TREE + "\n")
    return trees


# ---------------------------------------------------------------------
# edge_numbered_newick
# ---------------------------------------------------------------------
def test_edge_numbered_newick_labels_every_edge() -> None:
    """jplace requires each edge to carry a {N} tag so placements can reference it."""
    out = edge_numbered_newick(REF_TREE)
    # 4 tips + 2 internal + root = 7 edges in this topology
    assert out.count("{") == 7
    assert out.count("}") == 7
    assert out.endswith(";")


def test_edge_numbered_newick_numbers_are_unique_and_contiguous() -> None:
    """Duplicated or gapped edge numbers make the file unparseable."""
    nums = [int(n) for n in re.findall(r"\{(\d+)\}", edge_numbered_newick(REF_TREE))]
    assert sorted(nums) == list(range(len(nums)))


def test_edge_numbered_newick_preserves_tip_names_and_lengths() -> None:
    out = edge_numbered_newick(REF_TREE)
    for tip in ("A", "B", "C", "D"):
        assert tip in out
    assert "0.1" in out
    assert "0.35" in out


# ---------------------------------------------------------------------
# empty_jplace
# ---------------------------------------------------------------------
def test_empty_jplace_is_valid_json_with_the_real_tree() -> None:
    """An empty result must still carry the reference tree, not a placeholder.

    Downstream gappa commands read the tree to know what they are drawing; a
    fabricated stand-in would silently mislabel a figure.
    """
    doc = json.loads(empty_jplace(REF_TREE))
    assert doc["placements"] == []
    assert doc["version"] == 3
    assert "fields" in doc
    for tip in ("A", "B", "C", "D"):
        assert tip in doc["tree"]
    assert "{" in doc["tree"]  # edge-numbered


def test_empty_jplace_records_why_it_is_empty() -> None:
    """Zero placements and 'the stage never ran' are different claims."""
    doc = json.loads(empty_jplace(REF_TREE, reason="no queries for POL"))
    assert "no queries for POL" in json.dumps(doc["metadata"])


# ---------------------------------------------------------------------
# export_placement
# ---------------------------------------------------------------------
def test_export_copies_real_artifacts(tmp_path: Path) -> None:
    """The normal path: EPA-ng ran, so copy what it produced."""
    workdir = tmp_path / "wd" / "place_POL"
    workdir.mkdir(parents=True)
    (workdir / "epa_result.jplace").write_text('{"version":3,"placements":[1]}')
    (workdir / "labelled_tree.newick").write_text(REF_TREE)
    ref_dir = tmp_path / "ref"
    _write_ref_package(ref_dir)
    out = tmp_path / "out"

    jplace, newick = export_placement(
        workdir, ref_dir, "POL", out, stem="Toyus.ltr-flanked.POL"
    )

    assert jplace == out / "Toyus.ltr-flanked.POL.jplace"
    assert newick == out / "Toyus.ltr-flanked.POL.labelled.newick"
    assert json.loads(jplace.read_text())["placements"] == [1]
    assert newick.read_text().strip() == REF_TREE


def test_export_synthesises_an_empty_jplace_when_placement_did_not_run(
    tmp_path: Path,
) -> None:
    """No queries for this gene: the workdir does not even exist.

    A Snakemake rule declaring this output still needs the file, so the export
    writes a valid empty one rather than letting the DAG break.
    """
    ref_dir = tmp_path / "ref"
    _write_ref_package(ref_dir)
    out = tmp_path / "out"

    jplace, newick = export_placement(
        tmp_path / "never_created" / "place_POL",
        ref_dir,
        "POL",
        out,
        stem="Toyus.orphan.POL",
    )

    assert jplace.is_file()
    assert newick.is_file()
    doc = json.loads(jplace.read_text())
    assert doc["placements"] == []
    assert "A" in doc["tree"]  # the real reference tree, not a stand-in
    assert newick.read_text().strip() == REF_TREE


def test_export_prefers_the_raxml_optimised_tree(tmp_path: Path) -> None:
    """place() uses raxml.bestTree when present, so the export must match it.

    Exporting a different tree from the one placement actually ran on would make
    the published evidence disagree with the calls it produced.
    """
    ref_dir = tmp_path / "ref"
    trees = _write_ref_package(ref_dir)
    (trees / "POL.treefile").write_text("((A:9,B:9):9,(C:9,D:9):9);\n")  # stale
    out = tmp_path / "out"

    _, newick = export_placement(
        tmp_path / "absent" / "place_POL", ref_dir, "POL", out, stem="Toyus.x.POL"
    )
    assert "0.1" in newick.read_text()  # from raxml.bestTree, not the stale treefile


def test_export_raises_when_no_reference_tree_exists(tmp_path: Path) -> None:
    """Without any tree there is nothing honest to write; fail rather than fake."""
    out = tmp_path / "out"
    with pytest.raises(FileNotFoundError, match="POL"):
        export_placement(
            tmp_path / "wd", tmp_path / "empty_ref", "POL", out, stem="Toyus.x.POL"
        )


def test_export_creates_the_output_directory(tmp_path: Path) -> None:
    ref_dir = tmp_path / "ref"
    _write_ref_package(ref_dir)
    out = tmp_path / "deep" / "nested" / "out"
    jplace, _ = export_placement(
        tmp_path / "absent", ref_dir, "POL", out, stem="Toyus.x.POL"
    )
    assert jplace.is_file()
