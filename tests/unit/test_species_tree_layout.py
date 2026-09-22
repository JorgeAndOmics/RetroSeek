"""Tests for the early host species tree layout.

The tips it writes are joined, character for character, against the species names
the R plots show. If the two disagree, no row finds its tip and the tree silently
detaches. So the naming rule is pinned here against the R rule in style.R
(display_species): config value when present, else the stem with underscores as
spaces.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest
import species_tree_layout

pytest.importorskip("Bio")


def _config(tmp_path: Path, species: dict[str, str]) -> Path:
    lines = ["species:"] + [f"  {k}: '{v}'" for k, v in species.items()]
    p = tmp_path / "config.yaml"
    p.write_text("\n".join(lines) + "\n")
    return p


def _tips(out_dir: Path) -> list[dict[str, str]]:
    with (out_dir / "species.tree_tips.csv").open() as handle:
        return list(csv.DictReader(handle))


def test_display_name_prefers_the_config_value() -> None:
    assert (
        species_tree_layout.display_name("Homo_sapiens", {"Homo_sapiens": "Human"})
        == "Human"
    )


def test_display_name_falls_back_to_spaces_not_the_file_name() -> None:
    """Mirrors display_species() in style.R: never show a stem with underscores."""
    assert species_tree_layout.display_name("Myotis_myotis", {}) == "Myotis myotis"


def test_tips_carry_display_names_even_when_the_tree_uses_stems(tmp_path: Path) -> None:
    nwk = tmp_path / "t.nwk"
    nwk.write_text("((Homo_sapiens,Mus_musculus),Desmodus_rotundus);")
    cfg = _config(
        tmp_path,
        {
            "Homo_sapiens": "Homo sapiens",
            "Mus_musculus": "Mus musculus",
            "Desmodus_rotundus": "Desmodus rotundus",
        },
    )
    species_tree_layout.main(
        [
            "--config",
            str(cfg),
            "--genomes",
            "Homo_sapiens",
            "Mus_musculus",
            "Desmodus_rotundus",
            "--species-tree",
            str(nwk),
            "--out-dir",
            str(tmp_path / "out"),
        ]
    )
    assert sorted(t["tip"] for t in _tips(tmp_path / "out")) == [
        "Desmodus rotundus",
        "Homo sapiens",
        "Mus musculus",
    ]


def test_no_tree_configured_writes_header_only_files(tmp_path: Path) -> None:
    """The plots read header-only as "no tree" and fall back to config order."""
    cfg = _config(tmp_path, {"Homo_sapiens": "Homo sapiens"})
    species_tree_layout.main(
        [
            "--config",
            str(cfg),
            "--genomes",
            "Homo_sapiens",
            "--out-dir",
            str(tmp_path / "out"),
        ]
    )
    assert _tips(tmp_path / "out") == []
    assert (tmp_path / "out" / "species.tree_segments.csv").exists()
