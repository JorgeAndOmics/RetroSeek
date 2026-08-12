"""Unit tests for ``workflow/scripts/taxonomy/placement_figures.py``.

Turns the published placement evidence into per-(genome, tier) figures and
quality tables via gappa. Commands are built as argv lists and executed
separately, so these tests assert on the constructed command without needing
gappa installed - the same split ``run_ltr_retriever.py`` uses.

The load-bearing case is the empty jplace. ``gappa examine heat-tree`` does not
merely return nothing on a placement-free file: it aborts with "Invalid Color
Normalization with min >= max" and dumps core, because it tries to build a
colour scale from an empty mass range. Since a genome can legitimately place
nothing, the empty case is detected up front and an empty-state SVG is written
instead, mirroring the ``empty_plot()`` idiom the R plot generators already use.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from placement_figures import (
    count_placements,
    edpl_cmd,
    heat_tree_cmd,
    lwr_histogram_cmd,
    write_empty_state_svg,
)

REF_TREE = "((A:0.1{0},B:0.2{1}):0.3{2},(C:0.15{3},D:0.25{4}):0.35{5}){6};"


def _jplace(path: Path, n: int) -> Path:
    """A jplace with ``n`` placements on a fixed 4-tip tree."""
    path.write_text(
        json.dumps(
            {
                "version": 3,
                "fields": [
                    "edge_num",
                    "likelihood",
                    "like_weight_ratio",
                    "distal_length",
                    "pendant_length",
                ],
                "tree": REF_TREE,
                "placements": [
                    {"p": [[i % 6, -100.0, 1.0, 0.01, 0.02]], "n": [f"q{i}"]}
                    for i in range(n)
                ],
            }
        )
    )
    return path


# ---------------------------------------------------------------------
# count_placements
# ---------------------------------------------------------------------
def test_count_placements_counts_queries(tmp_path: Path) -> None:
    assert count_placements(_jplace(tmp_path / "a.jplace", 7)) == 7


def test_count_placements_of_an_empty_file_is_zero(tmp_path: Path) -> None:
    """This is the value that must gate every gappa call."""
    assert count_placements(_jplace(tmp_path / "e.jplace", 0)) == 0


def test_count_placements_of_a_missing_file_is_zero(tmp_path: Path) -> None:
    assert count_placements(tmp_path / "absent.jplace") == 0


# ---------------------------------------------------------------------
# command construction
# ---------------------------------------------------------------------
def test_heat_tree_cmd_requests_every_vector_format(tmp_path: Path) -> None:
    """SVG for editing, Newick for reuse, Nexus because FigTree opens it."""
    cmd = heat_tree_cmd(
        tmp_path / "x.jplace", tmp_path / "out", "Toyus.ltr-flanked.POL"
    )
    assert cmd[:3] == ["gappa", "examine", "heat-tree"]
    for flag in ("--write-svg-tree", "--write-newick-tree", "--write-nexus-tree"):
        assert flag in cmd
    assert "--allow-file-overwriting" in cmd


def test_heat_tree_cmd_carries_the_stem_as_file_prefix(tmp_path: Path) -> None:
    """gappa names files after the command; the prefix keeps runs distinguishable."""
    cmd = heat_tree_cmd(tmp_path / "x.jplace", tmp_path / "out", "Toyus.orphan.POL")
    assert "--file-prefix" in cmd
    assert cmd[cmd.index("--file-prefix") + 1] == "Toyus.orphan.POL."


def test_heat_tree_cmd_honours_mass_normalisation(tmp_path: Path) -> None:
    """`relative` makes genomes comparable; `absolute` preserves raw load."""
    cmd = heat_tree_cmd(
        tmp_path / "x.jplace", tmp_path / "o", "s", mass_norm="relative"
    )
    assert cmd[cmd.index("--mass-norm") + 1] == "relative"


def test_heat_tree_cmd_rejects_an_unknown_normalisation(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="mass_norm"):
        heat_tree_cmd(tmp_path / "x.jplace", tmp_path / "o", "s", mass_norm="banana")


def test_edpl_cmd_targets_the_right_subcommand(tmp_path: Path) -> None:
    """EDPL measures how spread out a query's placements are - an uncertainty
    axis independent of the confidence already in the catalog."""
    cmd = edpl_cmd(tmp_path / "x.jplace", tmp_path / "out", "Toyus.ltr-flanked.POL")
    assert cmd[:3] == ["gappa", "examine", "edpl"]
    assert "--file-prefix" in cmd


def test_lwr_histogram_cmd_targets_the_right_subcommand(tmp_path: Path) -> None:
    cmd = lwr_histogram_cmd(tmp_path / "x.jplace", tmp_path / "out", "s")
    assert cmd[:3] == ["gappa", "examine", "lwr-histogram"]


# ---------------------------------------------------------------------
# empty-state handling
# ---------------------------------------------------------------------
def test_write_empty_state_svg_produces_a_valid_standalone_svg(tmp_path: Path) -> None:
    """Snakemake declares this file, so it must exist even with nothing to draw."""
    out = tmp_path / "s.heat-tree.svg"
    write_empty_state_svg(out, "Toyus.orphan.POL", reason="no placements")
    text = out.read_text()
    assert text.lstrip().startswith("<svg") or "<svg" in text.split("\n")[0:3][0]
    assert "</svg>" in text


def test_empty_state_svg_says_which_sample_and_why(tmp_path: Path) -> None:
    """A blank figure with no explanation is indistinguishable from a broken one."""
    out = tmp_path / "s.svg"
    write_empty_state_svg(out, "Toyus.orphan.POL", reason="no placements for POL")
    text = out.read_text()
    assert "Toyus.orphan.POL" in text
    assert "no placements for POL" in text


def test_empty_state_svg_creates_parent_directories(tmp_path: Path) -> None:
    out = tmp_path / "deep" / "nested" / "s.svg"
    write_empty_state_svg(out, "x", reason="none")
    assert out.is_file()
