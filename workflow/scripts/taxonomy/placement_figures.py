"""Per-genome phylogenetic figures and quality tables from placement evidence.

Consumes the published `.jplace` files (see `taxonomy_placement.export_placement`)
and turns each into figures that answer questions the catalog cannot:

``heat-tree`` (written by gappa as ``{stem}.tree.svg``)
    The reference retroviral phylogeny with each branch coloured by the
    placement mass that landed on it. Reads as **where on the retroviral tree
    this genome's ERV load concentrates** - which lineages invaded it heavily,
    and which are represented by a handful of stragglers. Written as SVG (for
    editing), Newick (for reuse) and Nexus (FigTree opens it directly).

``edpl``
    Expected Distance between Placement Locations, per query. A locus whose
    placements are spread across distant branches is uncertain in a way its
    likelihood weight alone does not show, so this is an **independent
    uncertainty axis** from the ``confidence`` column already in the catalog. A
    locus with confidence 1.000 and a high EDPL is one whose certainty is an
    artifact of the LCA collapsing genuinely scattered placements.

``lwr-histogram``
    Distribution of likelihood weight ratios across the genome - a diagnostic
    for whether the reference tree resolves this material at all.

Deliberately no ``gappa examine graft``: it draws one pendant edge per query and
is unreadable past a few hundred (Mus musculus places 5,690 in the LTR-flanked
tier alone). Heat-trees accumulate mass onto branches and stay legible at any
count.

The empty case matters
----------------------
``gappa examine heat-tree`` does not return quietly on a placement-free file. It
aborts with ``Invalid Color Normalization with min >= max`` and dumps core,
because it builds a colour scale from an empty mass range. A genome can
legitimately place nothing, so every gappa call is gated on the placement count
and an empty-state SVG is written instead - the same idea as ``empty_plot()`` in
the R plot generators.

CLI
---
::

    python placement_figures.py \\
        --jplace <path> --out-dir <dir> --stem <genome>.<tier>.<gene> \\
        [--mass-norm absolute|relative] [--skip-edpl] [--skip-lwr]
"""

from __future__ import annotations

import argparse
import json
import logging
import shutil
from pathlib import Path

from external import run_tool
from log import OK, job_logging, run_main

logger = logging.getLogger(__name__)

GAPPA = "gappa"  # resolved from PATH (the RetroSeek conda env)
VALID_MASS_NORM = ("absolute", "relative")

# The house colours this script needs (docs/visual_style.md). plot2sort/style.R
# is the single source of truth; these mirror it, and
# test_house_colours_mirror_style_r fails if the two drift apart.
PAPER, INK, INK_SOFT = "#FFFFFF", "#222222", "#5A5A5A"
# Mid grey to indigo: a branch with no placement mass must still be visible, or
# the tree loses its shape, so the ramp starts at the house mid grey rather than
# the near-white light end of the R figures' ramp.
RAMP = ("#BBBBBB", "#332288")


def count_placements(jplace: Path) -> int:
    """Number of placed queries in a jplace file; 0 if absent or unreadable.

    This is the gate on every gappa call. Treating a missing or malformed file
    as zero is deliberate: the caller then writes an empty-state figure rather
    than crashing, and the reason is logged.
    """
    if not jplace.is_file():
        return 0
    try:
        doc = json.loads(jplace.read_text())
    except (json.JSONDecodeError, OSError):
        logger.warning("could not parse %s; treating as empty", jplace)
        return 0
    placements = doc.get("placements", [])
    return len(placements) if isinstance(placements, list) else 0


def _base_cmd(sub: str, jplace: Path, out_dir: Path, stem: str) -> list[str]:
    """Shared gappa argv: subcommand, input, output dir, and the naming prefix.

    gappa names outputs after the command it ran, so without a prefix every
    genome would overwrite the last one's files in a shared directory.
    """
    return [
        GAPPA,
        "examine",
        sub,
        "--jplace-path",
        str(jplace),
        "--out-dir",
        str(out_dir),
        "--file-prefix",
        f"{stem}.",
        "--allow-file-overwriting",
    ]


def heat_tree_cmd(
    jplace: Path, out_dir: Path, stem: str, mass_norm: str = "absolute"
) -> list[str]:
    """Build the heat-tree argv.

    ``mass_norm`` mirrors gappa's own option: ``absolute`` keeps raw placement
    mass, so a genome with more ERVs looks hotter; ``relative`` normalises
    within each sample, which is what you want when comparing genomes of very
    different ERV load.
    """
    if mass_norm not in VALID_MASS_NORM:
        raise ValueError(
            f"unknown mass_norm {mass_norm!r}; expected one of {VALID_MASS_NORM}"
        )
    return [
        *_base_cmd("heat-tree", jplace, out_dir, stem),
        "--mass-norm",
        mass_norm,
        # House colours: no mass in mid grey, the heaviest branches in indigo.
        "--color-list",
        ",".join(RAMP),
        "--write-svg-tree",
        "--write-newick-tree",
        "--write-nexus-tree",
    ]


def edpl_cmd(jplace: Path, out_dir: Path, stem: str) -> list[str]:
    """Build the EDPL argv (per-query placement spread + its histogram)."""
    return _base_cmd("edpl", jplace, out_dir, stem)


def lwr_histogram_cmd(jplace: Path, out_dir: Path, stem: str) -> list[str]:
    """Build the likelihood-weight-ratio histogram argv."""
    return _base_cmd("lwr-histogram", jplace, out_dir, stem)


def write_empty_state_svg(path: Path, stem: str, reason: str) -> None:
    """Write a standalone SVG standing in for a figure that has no data.

    Names the sample and the reason: a blank figure with no explanation is
    indistinguishable from a broken one, and the two need different responses.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" width="640" height="140" '
        'viewBox="0 0 640 140">\n'
        f'  <rect width="640" height="140" fill="{PAPER}"/>\n'
        '  <text x="320" y="62" text-anchor="middle" '
        'font-family="IBM Plex Sans, sans-serif" '
        f'font-size="15" font-weight="bold" fill="{INK}">{stem}</text>\n'
        '  <text x="320" y="90" text-anchor="middle" '
        'font-family="IBM Plex Sans, sans-serif" '
        f'font-size="13" fill="{INK_SOFT}">{reason}</text>\n'
        "</svg>\n"
    )


def route_heat_tree_outputs(plot_dir: Path, tree_dir: Path, stem: str) -> None:
    """Move the heat-tree's Newick and Nexus out of the figure directory.

    gappa writes an SVG, a Newick and a Nexus from one command into a single
    --out-dir. Only the SVG is a figure; the other two are tree artifacts and
    belong beside the .jplace they came from. Missing files are ignored - the
    empty-placement path writes a placeholder SVG and nothing else.
    """
    tree_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("tree.newick", "tree.nexus"):
        src = plot_dir / f"{stem}.{ext}"
        if src.is_file():
            shutil.move(str(src), str(tree_dir / src.name))


def main(argv: list[str] | None = None) -> None:
    """Entry point."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--jplace", type=Path, required=True)
    parser.add_argument(
        "--out-dir", type=Path, required=True, help="figures (SVG) land here"
    )
    parser.add_argument(
        "--table-dir",
        type=Path,
        default=None,
        help="EDPL and LWR CSVs; defaults to --out-dir when unset",
    )
    parser.add_argument(
        "--tree-dir",
        type=Path,
        default=None,
        help="heat-tree Newick/Nexus; defaults to --out-dir when unset",
    )
    parser.add_argument(
        "--stem", required=True, help="{genome}.{tier}.{gene}; names every output"
    )
    parser.add_argument("--mass-norm", choices=VALID_MASS_NORM, default="absolute")
    parser.add_argument("--skip-edpl", action="store_true")
    parser.add_argument("--skip-lwr", action="store_true")
    parser.add_argument("--log", type=Path, help="job log (the Snakemake log: path)")
    args = parser.parse_args(argv)
    job_logging(args.log, "placement_figures")
    table_dir = args.table_dir or args.out_dir
    tree_dir = args.tree_dir or args.out_dir
    for d in (args.out_dir, table_dir, tree_dir):
        d.mkdir(parents=True, exist_ok=True)

    n = count_placements(args.jplace)
    if n == 0:
        # Every gappa command here divides by a mass range that does not exist
        # when nothing was placed; heat-tree in particular core-dumps.
        logger.warning(
            "%s has no placements; writing empty-state figures instead", args.stem
        )
        write_empty_state_svg(
            args.out_dir / f"{args.stem}.tree.svg",
            args.stem,
            reason="no queries were placed on the reference tree",
        )
        return

    logger.info("%s: %d placed queries", args.stem, n)
    run_tool(heat_tree_cmd(args.jplace, args.out_dir, args.stem, args.mass_norm))
    route_heat_tree_outputs(args.out_dir, tree_dir, args.stem)
    if not args.skip_edpl:
        run_tool(edpl_cmd(args.jplace, table_dir, args.stem))
    if not args.skip_lwr:
        run_tool(lwr_histogram_cmd(args.jplace, table_dir, args.stem))
    logger.log(OK, "%s placed queries drawn", f"{n:,}")


if __name__ == "__main__":
    run_main(main)
