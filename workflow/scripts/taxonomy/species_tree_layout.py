"""Lay out the host species tree from the config alone, early in the pipeline.

Every figure that compares genomes puts them on rows in one canonical order, the
host tree's when one is configured (docs/visual_style.md). That order has to be
available to every stage, including the ranges-stage plots that run long before
classification. So the species tree is laid out here from nothing but the user's
Newick and the config, rather than inside the taxonomy stage from classified-locus
tables as it used to be.

The tips are the configured genomes under their DISPLAY names, exactly as the R
side's display_species() writes them: the config `species:` value when there is
one, otherwise the file stem with underscores as spaces. The two must agree
character for character or no row would find its tip. A tree labelled with stems
still matches, because the stems are passed as aliases.

With no tree configured the files are written header-only, which the plotting code
reads as "no tree": species then fall back to config order, and the DAG is stable
either way.

Outputs, in `--out-dir`: species.tree_tips.csv (tip, x, y) and
species.tree_segments.csv (x, y, xend, yend).
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import yaml
from tree_layout import build_alias_index, from_newick, layout, write

logger = logging.getLogger(__name__)


def display_name(stem: str, species_map: dict[str, str]) -> str:
    """The name a plot shows for a genome; mirrors display_species() in style.R."""
    value = species_map.get(stem)
    return str(value) if value else stem.replace("_", " ")


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--genomes",
        nargs="+",
        required=True,
        help="genome stems in the study, in config order",
    )
    parser.add_argument("--species-tree", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)

    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8")) or {}
    species_map = {str(k): str(v) for k, v in (cfg.get("species") or {}).items()}
    tips = [display_name(stem, species_map) for stem in args.genomes]

    if args.species_tree and str(args.species_tree) and args.species_tree.exists():
        tree = from_newick(
            args.species_tree, tips, aliases=build_alias_index(species_map)
        )
        # Real branch lengths are kept; a bare cladogram is squared off, which
        # reads better beside rows.
        clades = tree.find_clades()  # type: ignore[no-untyped-call]
        has_lengths = any(c.branch_length not in (None, 0) for c in clades)
        write(args.out_dir, "species", *layout(tree, align_tips=not has_lengths))
        n_tips = tree.count_terminals()  # type: ignore[no-untyped-call]
        logger.info("species tree: %d tips laid out", n_tips)
    else:
        logger.info("species tree: none configured, species will follow config order")
        write(args.out_dir, "species", [], [])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
