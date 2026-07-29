"""
tree_layout.py
==============
Turns the two trees RetroSeek can draw into flat COORDINATE tables that the R
plot generators render with `geom_segment` — no R tree library required.

Two sources, deliberately different in kind:

* **taxon tree** — built from the reference's ``taxonomy.tsv`` (``name, parent,
  rank``). This is the ADR-008 classification axis rendered as a cladogram: it
  covers every axis taxon by construction, needs no representative-sequence
  choice, and is rank-agnostic. It is a *classification*, not a phylogeny, so it
  carries no branch lengths (ADR-011).
* **species tree** — a user-supplied Newick pinned by ``input.species_tree``.
  Reproducible (no network at run time) and free to carry real divergence times
  if the user exports a dated tree.

Parsing/pruning uses ``Bio.Phylo`` from biopython, which the pipeline already
pins for the tBLASTn cache — so trees cost ZERO new dependencies.

Output per tree, written next to the classification tables:
  ``<name>.tree_segments.csv``  x, y, xend, yend   (the drawn lines)
  ``<name>.tree_tips.csv``      tip, x, y          (leaf label anchors)

Determinism: tips are ladderized and every sibling set is ordered by a stable
name key, so the same inputs always produce byte-identical coordinates.
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq
import yaml
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

logger = logging.getLogger(__name__)

# Written even when a tree is unavailable, so the Snakemake DAG stays stable and
# the R side renders an explanatory placeholder instead of failing.
SEGMENT_HEADER = ("x", "y", "xend", "yend")
TIP_HEADER = ("tip", "x", "y")


# --------------------------------------------------------------------- loading
def load_hierarchy(taxonomy_tsv: Path) -> tuple[dict[str, str | None], dict[str, str]]:
    """Return (child -> parent, name -> rank) from a ``name,parent,rank`` TSV.

    Same file `taxonomy_lca.load_taxonomy` consumes; read independently here so
    the layout step has no import-time coupling to the classifier.
    """
    parent: dict[str, str | None] = {}
    rank: dict[str, str] = {}
    with taxonomy_tsv.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            parent[row["name"]] = row["parent"] or None
            rank[row["name"]] = row["rank"]
    return parent, rank


def _ancestry(node: str, parent: dict[str, str | None]) -> list[str]:
    """``node`` then every ancestor up to the root (cycle-safe)."""
    chain: list[str] = []
    seen: set[str] = set()
    cur: str | None = node
    while cur is not None and cur not in seen:
        chain.append(cur)
        seen.add(cur)
        cur = parent.get(cur)
    return chain


def from_taxonomy(taxonomy_tsv: Path, tips: list[str]) -> Tree | None:
    """Build a cladogram of ``tips`` from the taxonomy hierarchy.

    Only the observed taxa and the ancestors connecting them are kept, so the
    tree shows the study's actual axis rather than the whole reference.
    Returns None when nothing can be placed.
    """
    parent, _rank = load_hierarchy(taxonomy_tsv)
    known = [t for t in tips if t in parent]
    missing = sorted(set(tips) - set(known))
    if missing:
        logger.warning(
            "taxon tree: %d call(s) absent from taxonomy.tsv, omitted: %s",
            len(missing),
            ", ".join(missing[:10]),
        )
    if not known:
        return None

    # Keep every node on a path from an observed tip to its root.
    keep: set[str] = set()
    for t in known:
        keep.update(_ancestry(t, parent))

    children: dict[str | None, list[str]] = {}
    for node in sorted(keep):  # sorted => deterministic
        children.setdefault(
            parent.get(node) if parent.get(node) in keep else None, []
        ).append(node)

    tip_set = set(known)

    def build(name: str) -> Clade:
        kids = [k for k in children.get(name, []) if k != name]
        if not kids:
            return Clade(name=name, branch_length=1.0)
        sub = [build(k) for k in kids]
        # Mixed-rank axis (ADR-008): a taxon can be BOTH an observed call and an
        # ancestor of other calls - e.g. loci resolved only to `Retroviridae`
        # alongside loci resolved to `Gammaretrovirus`. It then needs its own tip
        # (so its bar has a row) *and* its internal node (so its descendants hang
        # off it). The internal node stays unnamed to avoid a duplicate label.
        if name in tip_set:
            sub.insert(0, Clade(name=name, branch_length=1.0))
        return Clade(name=None, branch_length=1.0, clades=sub)

    roots = children.get(None, [])
    if not roots:
        return None
    top = (
        build(roots[0])
        if len(roots) == 1
        else Clade(name=None, branch_length=1.0, clades=[build(r) for r in roots])
    )
    return Tree(root=top, rooted=True)


def _normalize(label: str) -> str:
    """Fold a tip label / species name to a comparison key.

    Newick exports differ from config display names in separators and case
    (``Mus_musculus`` vs ``Mus musculus``), so match on a normalized key.
    """
    return label.replace("_", " ").strip().lower()


def from_newick(newick: Path, tips: list[str]) -> Tree:
    """Read a user Newick and prune it to ``tips`` (matched name-insensitively).

    Reports unmatched names in BOTH directions: a species missing from the tree
    silently vanishes from the figure otherwise, and that is exactly the kind of
    silent data loss this pipeline fails loudly on elsewhere.
    """
    tree = Phylo.read(str(newick), "newick")
    want = {_normalize(t): t for t in tips}
    matched: dict[str, str] = {}  # tip label -> display name
    for leaf in tree.get_terminals():
        key = _normalize(leaf.name or "")
        if key in want:
            matched[leaf.name] = want[key]

    unmatched_species = sorted(set(tips) - set(matched.values()))
    if unmatched_species:
        logger.warning(
            "species tree: %d species not found among the tree tips: %s",
            len(unmatched_species),
            ", ".join(unmatched_species[:10]),
        )
    if not matched:
        raise SystemExit(
            f"species tree {newick} shares no tip with the study's species. "
            "Check that tip labels match the config `species:` display names."
        )

    for leaf in list(tree.get_terminals()):
        if leaf.name not in matched:
            tree.prune(leaf)
    for leaf in tree.get_terminals():
        leaf.name = matched.get(leaf.name, leaf.name)
    result: Tree = tree
    return result


# ---------------------------------------------------------------------- layout
def layout(
    tree: Tree, align_tips: bool
) -> tuple[list[tuple[float, ...]], list[tuple[Any, ...]]]:
    """Rectangular tree coordinates.

    x = cumulative branch length from the root (unit edges when the tree has no
    lengths); y = ladderized tip order; an internal node sits at the mean y of
    its children. ``align_tips`` squares the leaves off at the deepest x, which
    is what a cladogram wants so labels line up with the bars beside them.
    """
    tree.ladderize()
    for i, t in enumerate(tree.get_terminals(), 1):
        t.y = float(i)

    def set_x(clade: Clade, x0: float) -> None:
        bl = clade.branch_length
        bl = 1.0 if bl is None else float(bl)
        clade.x = x0 + (0.0 if clade is tree.root else bl)
        for child in clade.clades:
            set_x(child, clade.x)

    tree.root.x = 0.0
    for child in tree.root.clades:
        set_x(child, 0.0)

    def set_y(clade: Clade) -> float:
        if clade.is_terminal():
            return float(clade.y)
        ys = [set_y(c) for c in clade.clades]
        clade.y = sum(ys) / len(ys)
        return float(clade.y)

    set_y(tree.root)

    max_x = max(c.x for c in tree.find_clades())
    if align_tips:
        for t in tree.get_terminals():
            t.x = max_x

    segments: list[tuple[float, ...]] = []
    for clade in tree.find_clades():
        if clade.is_terminal():
            continue
        ys = [c.y for c in clade.clades]
        segments.append((clade.x, min(ys), clade.x, max(ys)))  # vertical spine
        segments.extend(  # horizontal arms
            (clade.x, child.y, child.x, child.y) for child in clade.clades
        )
    tips = [(t.name, t.x, t.y) for t in tree.get_terminals()]
    return segments, tips


def write(
    out_dir: Path,
    name: str,
    segments: list[tuple[float, ...]],
    tips: list[tuple[Any, ...]],
) -> None:
    """Write the two CSVs for one tree (headers always, rows possibly none)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    with (out_dir / f"{name}.tree_segments.csv").open(
        "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.writer(fh)
        w.writerow(SEGMENT_HEADER)
        w.writerows(segments)
    with (out_dir / f"{name}.tree_tips.csv").open(
        "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.writer(fh)
        w.writerow(TIP_HEADER)
        w.writerows(tips)
    logger.info("%s tree: %d tips, %d segments", name, len(tips), len(segments))


# ------------------------------------------------------------------------- CLI
def _observed(parquet_dir: Path, column: str) -> list[str]:
    """Distinct non-empty values of ``column`` across the loci/orphan tables.

    Reads the per-genome tables rather than catalog.csv on purpose: the catalog
    is an OUTPUT of the plot stage, so depending on it here would create a cycle.
    """
    seen: set[str] = set()
    for pattern in ("*.loci.parquet", "*.orphans.parquet"):
        for path in sorted(parquet_dir.glob(pattern)):
            table = pq.read_table(path)
            if column == "species":
                # `species` is not a table column - it is the filename stem.
                seen.add(path.name.split(".")[0])
                continue
            if column not in table.column_names:
                continue
            seen.update(str(v) for v in table.column(column).to_pylist() if v)
    return sorted(seen)


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--parquet-dir",
        type=Path,
        required=True,
        help="dir holding <genome>.loci.parquet / .orphans.parquet",
    )
    p.add_argument(
        "--taxonomy-tsv",
        type=Path,
        required=True,
        help="reference taxonomy.tsv (name, parent, rank)",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="destination for the *.tree_segments.csv / *.tree_tips.csv",
    )
    p.add_argument(
        "--species-tree",
        type=Path,
        default=None,
        help="optional user-supplied Newick of the host species",
    )
    p.add_argument(
        "--config",
        type=Path,
        default=None,
        help="pipeline config YAML; its `species:` map canonicalizes genome "
        "stems to the display names the plots use, so tree tips match the bars",
    )
    args = p.parse_args()

    # --- taxon tree (always available: taxonomy.tsv ships with the reference)
    taxa = _observed(args.parquet_dir, "taxon_call")
    tree = from_taxonomy(args.taxonomy_tsv, taxa) if taxa else None
    if tree is None:
        logger.warning("taxon tree: nothing to draw (no resolvable taxon calls)")
        write(args.out_dir, "taxon", [], [])
    else:
        write(args.out_dir, "taxon", *layout(tree, align_tips=True))

    # --- species tree (only when the user pinned one)
    # The plot generators relabel genome stems to the config `species:` display
    # names before plotting, so the tips must carry those same names or nothing
    # would match (`Homo_sapiens` tip vs `Homo sapiens` bar).
    species = _observed(args.parquet_dir, "species")
    if args.config and args.config.exists():
        cfg = yaml.safe_load(args.config.read_text(encoding="utf-8")) or {}
        mapping = cfg.get("species") or {}
        species = sorted({str(mapping.get(s, s)) for s in species})
    if args.species_tree and str(args.species_tree) and args.species_tree.exists():
        stree = from_newick(args.species_tree, species)
        # A user tree may carry real branch lengths; keep them (align_tips=False)
        # unless it is a bare cladogram, where squared-off tips read better.
        has_len = any(c.branch_length not in (None, 0) for c in stree.find_clades())
        write(args.out_dir, "species", *layout(stree, align_tips=not has_len))
    else:
        logger.info("species tree: none configured (input.species_tree unset)")
        write(args.out_dir, "species", [], [])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
