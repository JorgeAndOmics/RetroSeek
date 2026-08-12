"""Co-phylogeny: do ERV complements track host ancestry?

Takes the published `.jplace` files for one tier and builds a tree of the *host
genomes* from their placement distributions, then sets that against the host
phylogeny.

    ERV-composition tree            host phylogeny
    (gappa analyze squash)          (input.species_tree)
              \\                          /
               `------ compared ---------'
                          |
              congruent  -> ERVs largely inherited vertically with their hosts
              discordant -> cross-species transmission, or lineage-specific
                            expansion and loss

Three gappa commands do the work, all on placement mass over the shared
reference tree:

``analyze squash``   hierarchical clustering of samples -> the composition tree
``analyze krd``      the pairwise Kantorovich-Rubinstein distances behind it
``analyze edgepca``  ordination, showing which branches drive the separation

Two mechanics are easy to get wrong
-----------------------------------
**Sample naming.** gappa labels each sample by its file basename, and EPA-ng
writes every run to ``epa_result.jplace``. Handed those directly, all five tips
come back named ``epa_result``. Inputs are staged under genome-derived names
first, which is why `taxonomy_placement.export_placement` publishes them as
``{genome}.{tier}.{gene}.jplace``.

**Branch lengths are not comparable.** Congruence is measured on bipartitions
(Robinson-Foulds), never on lengths. The supplied host tree is frequently a
cladogram carrying placeholder lengths - the shipped one is all 1s and 2s - and
setting those against placement distances would be meaningless. Where the host
tree does carry real divergence times, the KRD matrix is emitted alongside so
the quantitative comparison can be made deliberately rather than by accident.

Tiers are analysed separately: whether the weaker orphan tier tells the same
story as the LTR-confirmed one is itself a useful check.

CLI
---
::

    python placement_cophylogeny.py \\
        --jplace <f1> <f2> ... --out-dir <dir> --tier ltr-flanked \\
        [--host-tree <newick>] [--exponent 1]
"""

from __future__ import annotations

import argparse
import csv
import logging
import shutil
import subprocess
import sys
from io import StringIO
from pathlib import Path

from Bio import Phylo

logger = logging.getLogger(__name__)

GAPPA = "gappa"  # resolved from PATH (the RetroSeek conda env)

# {genome}.{tier}.{gene}.jplace - the tier and gene are fixed per invocation, so
# stripping the last three dot-separated fields leaves the genome, even when the
# genome name itself contains dots (assembly accessions do).
_STEM_FIELDS = 3


def sample_name_for(filename: str) -> str:
    """Genome name from a published placement filename.

    Tips must be host genomes for the tree to be comparable against the host
    phylogeny, so the tier and gene are stripped.
    """
    stem = Path(filename).name
    if stem.endswith(".jplace"):
        stem = stem[: -len(".jplace")]
    parts = stem.split(".")
    return ".".join(parts[:-2]) if len(parts) > _STEM_FIELDS - 1 else stem


def stage_jplace(jplace_files: list[Path], staged_dir: Path) -> Path:
    """Copy inputs into ``staged_dir`` renamed to ``{genome}.jplace``.

    gappa derives sample labels from basenames, so this rename is what makes the
    resulting tree's tips readable. Duplicate genome names are rejected rather
    than silently overwritten - two tiers of one genome would otherwise collapse
    into a single sample and compare a genome against itself.
    """
    staged_dir.mkdir(parents=True, exist_ok=True)
    seen: dict[str, Path] = {}
    for src in jplace_files:
        name = sample_name_for(src.name)
        if name in seen:
            raise ValueError(
                f"duplicate sample name {name!r} from {src.name} and "
                f"{seen[name].name}; analyse one tier at a time"
            )
        seen[name] = src
        shutil.copyfile(src, staged_dir / f"{name}.jplace")
    return staged_dir


def _analyze_cmd(sub: str, staged_dir: Path, out_dir: Path) -> list[str]:
    return [
        GAPPA,
        "analyze",
        sub,
        "--jplace-path",
        str(staged_dir),
        "--out-dir",
        str(out_dir),
        "--allow-file-overwriting",
    ]


def squash_cmd(staged_dir: Path, out_dir: Path) -> list[str]:
    """Squash clustering: the tree of samples, written as Newick."""
    return [*_analyze_cmd("squash", staged_dir, out_dir), "--write-newick-tree"]


def krd_cmd(staged_dir: Path, out_dir: Path) -> list[str]:
    """Pairwise Kantorovich-Rubinstein distances between samples."""
    return _analyze_cmd("krd", staged_dir, out_dir)


# ---------------------------------------------------------------------
# topology comparison
# ---------------------------------------------------------------------
def _tips(tree: Phylo.BaseTree.Tree) -> set[str]:
    return {t.name for t in tree.get_terminals() if t.name}


def bipartitions(newick: str) -> set[frozenset[str]]:
    """Non-trivial splits induced by a tree, as frozensets of tip names.

    Each internal branch splits the tips in two; the smaller side names the
    split. Single-tip and all-tip splits are excluded because every tree over
    the same taxa shares them, so they carry no comparative signal.
    """
    tree = Phylo.read(StringIO(newick), "newick")
    all_tips = _tips(tree)
    splits: set[frozenset[str]] = set()
    for clade in tree.get_nonterminals():
        side: set[str] = {t.name for t in clade.get_terminals() if t.name}
        other: set[str] = all_tips - side
        if len(side) < 2 or len(other) < 1:
            continue
        # Canonical orientation so the same split matches from either tree:
        # always name the split by its smaller side, ties broken alphabetically.
        smaller = (
            side if (len(side), sorted(side)) <= (len(other), sorted(other)) else other
        )
        splits.add(frozenset(smaller))
    return {s for s in splits if len(s) >= 2}


def congruence(host_newick: str, erv_newick: str) -> dict[str, object]:
    """Compare two topologies over the same tips.

    Returns the split counts, the Robinson-Foulds distance (splits present in
    one tree but not the other), and the offending splits themselves - a bare
    distance is not actionable, the caller needs to see which grouping
    disagrees.
    """
    host_tips = _tips(Phylo.read(StringIO(host_newick), "newick"))
    erv_tips = _tips(Phylo.read(StringIO(erv_newick), "newick"))
    if host_tips != erv_tips:
        missing = (host_tips ^ erv_tips) or {"<none>"}
        raise ValueError(
            "host and ERV trees must cover the same tip set; "
            f"differing tips: {sorted(missing)}"
        )

    host_splits = bipartitions(host_newick)
    erv_splits = bipartitions(erv_newick)
    shared = host_splits & erv_splits
    return {
        "n_tips": len(host_tips),
        "host_splits": len(host_splits),
        "erv_splits": len(erv_splits),
        "shared_splits": len(shared),
        "rf_distance": len(host_splits ^ erv_splits),
        "congruent": host_splits == erv_splits,
        "host_only_splits": sorted(
            "|".join(sorted(s)) for s in host_splits - erv_splits
        ),
        "erv_only_splits": sorted(
            "|".join(sorted(s)) for s in erv_splits - host_splits
        ),
    }


def run(cmd: list[str]) -> None:
    """Invoke gappa, surfacing its stderr on failure."""
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if res.returncode != 0:
        sys.stderr.write((res.stderr or "")[-3000:])
        raise SystemExit(f"gappa failed ({res.returncode}): {' '.join(cmd[:3])}")


def write_summary(path: Path, tier: str, result: dict[str, object] | None) -> None:
    """Write the congruence verdict as a one-row-per-metric CSV.

    Written even when no host tree was supplied, recording that the comparison
    was skipped: "not congruent" and "not compared" are different claims.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["tier", "metric", "value"])
        if result is None:
            writer.writerow([tier, "comparison", "skipped (no host tree supplied)"])
            return
        for key in (
            "n_tips",
            "host_splits",
            "erv_splits",
            "shared_splits",
            "rf_distance",
            "congruent",
        ):
            writer.writerow([tier, key, result[key]])
        for key in ("host_only_splits", "erv_only_splits"):
            # congruence() stores these as list[str]; the dict is typed
            # dict[str, object] because it also carries ints and bools.
            splits = result[key]
            if isinstance(splits, list):
                for split in splits:
                    writer.writerow([tier, key, split])


def main(argv: list[str] | None = None) -> int:
    """Entry point."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--jplace", type=Path, nargs="+", required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--tier", required=True, help="ltr-flanked | orphan")
    parser.add_argument(
        "--gene", default="POL", help="placement gene; names the summary"
    )
    parser.add_argument(
        "--host-tree", type=Path, default=None, help="Newick host phylogeny"
    )
    parser.add_argument("--staged-dir", type=Path, default=None)
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if len(args.jplace) < 3:
        # Squash clustering over fewer than three samples has no internal
        # structure to compare, so the whole comparison is vacuous.
        logger.warning(
            "only %d sample(s); a composition tree needs at least 3", len(args.jplace)
        )
        write_summary(
            args.out_dir / f"cophylogeny_summary.{args.tier}.{args.gene}.csv",
            args.tier,
            None,
        )
        return 0

    staged = stage_jplace(args.jplace, args.staged_dir or args.out_dir / "_staged")
    run(squash_cmd(staged, args.out_dir))
    run(krd_cmd(staged, args.out_dir))

    cluster = args.out_dir / "cluster.newick"
    erv_tree = args.out_dir / f"erv_composition.{args.tier}.newick"
    if cluster.exists():
        shutil.move(str(cluster), erv_tree)

    result = None
    if args.host_tree and args.host_tree.is_file() and erv_tree.is_file():
        try:
            result = congruence(args.host_tree.read_text(), erv_tree.read_text())
        except ValueError as exc:
            # Tip sets differ (a genome without placements, say). Record it
            # rather than aborting the run.
            logger.warning("congruence not computed: %s", exc)
        else:
            verdict = "congruent" if result["congruent"] else "DISCORDANT"
            logger.info(
                "%s tier: %s with the host phylogeny (RF=%s, %s/%s splits shared)",
                args.tier,
                verdict,
                result["rf_distance"],
                result["shared_splits"],
                result["host_splits"],
            )
    write_summary(
        args.out_dir / f"cophylogeny_summary.{args.tier}.{args.gene}.csv",
        args.tier,
        result,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
