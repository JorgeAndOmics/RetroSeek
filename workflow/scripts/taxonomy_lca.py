"""
Taxonomy LCA for ERV locus classification
==========================================

Resolves the taxonomic identity of an ERV locus from the *set* of retroviral
genera supported by its evidence, using a lowest-common-ancestor (LCA) walk over
a small curated retroviral taxonomy.

Motivation
----------
RetroSeek detection produces, per locus, a cloud of cross-genus probe hits
(retroviral ``pol``/RT is conserved across genera, so one genomic ERV is
homologous to many reference viruses at once). The legacy taxon assignment either
picks the single highest-bitscore probe (``best`` — an uncalibrated argmax) or
dumps the whole genus set (``concatenate`` — uninterpretable). Neither classifies.

This module replaces that with LCA backoff (the MEGAN / detectEVE paradigm,
expressed over a taxonomy rather than a tree):

* genus set collapses to a single genus  -> assign that **genus** (confident);
* genus set spans several genera          -> assign their **LCA** (subfamily /
  family — an honest higher rank);
* nothing retroviral supported            -> ``unclassified``.

It is intentionally probe-agnostic: it operates on whatever genus labels the
evidence carries, never on a fixed marker.

This is the *unweighted* baseline (presence/absence of genus support). The
weighted variant (blastx bitscore-weighted LCA against an independent reference)
builds on the same taxonomy and ``lca``/``rank_of`` primitives — see
``docs/taxonomy_classification/``.
"""

from __future__ import annotations

import csv
import re
import sys
from collections import Counter
from pathlib import Path

# =============================================================================
# 1. Retroviral taxonomy (parent pointers)
# =============================================================================
# The OPERATIONAL hierarchy is data-derived: callers run ``load_taxonomy()`` to
# replace the maps below with one built from NCBI Taxonomy
# (``taxonomy_build_hierarchy.py`` -> ``taxonomy.tsv``), so the hierarchy is not
# hard-coded and self-updates from the authoritative source. The literals here are
# only a FALLBACK DEFAULT (unit tests / when no taxonomy.tsv is present).
# Parent map: child -> parent. Root has parent None.
RETRO_PARENT: dict[str, str | None] = {
    "Retroviridae": None,
    "Orthoretrovirinae": "Retroviridae",
    "Spumaretrovirinae": "Retroviridae",
    "Alpharetrovirus": "Orthoretrovirinae",
    "Betaretrovirus": "Orthoretrovirinae",
    "Gammaretrovirus": "Orthoretrovirinae",
    "Deltaretrovirus": "Orthoretrovirinae",
    "Epsilonretrovirus": "Orthoretrovirinae",
    "Lentivirus": "Orthoretrovirinae",
    "Spumaretrovirus": "Spumaretrovirinae",
}

RANK_OF: dict[str, str] = {
    "Retroviridae": "family",
    "Orthoretrovirinae": "subfamily",
    "Spumaretrovirinae": "subfamily",
    "Alpharetrovirus": "genus",
    "Betaretrovirus": "genus",
    "Gammaretrovirus": "genus",
    "Deltaretrovirus": "genus",
    "Epsilonretrovirus": "genus",
    "Lentivirus": "genus",
    "Spumaretrovirus": "genus",
}

# ERV class system (Jern/Blomberg) — a curated biological grouping (NOT an NCBI rank,
# so it can't be derived from taxonomy). Operationally loaded from an overridable data
# file via ``load_erv_class()``; the literal below is only a fallback default.
ERV_CLASS: dict[str, str] = {
    "Gammaretrovirus": "Class I",
    "Epsilonretrovirus": "Class I",
    "Alpharetrovirus": "Class II",
    "Betaretrovirus": "Class II",
    "Deltaretrovirus": "Class II",
    "Lentivirus": "Class II",
    "Spumaretrovirus": "Class III",
}

UNCLASSIFIED = "unclassified"


def load_taxonomy(path: str | Path) -> dict[str, str | None]:
    """
    Replace the hierarchy with a data-derived one (``taxonomy.tsv``).

    Columns: ``name, parent, rank`` (root's parent empty). Built by
    ``taxonomy_build_hierarchy.py`` from NCBI Taxonomy. Reassigns the module-level
    ``RETRO_PARENT``/``RANK_OF`` so ``ancestors``/``lca``/``rank_of`` use it. The
    literal maps above are only a fallback when this is not called.
    """
    global RETRO_PARENT, RANK_OF
    parent: dict[str, str | None] = {}
    rank: dict[str, str] = {}
    with Path(path).open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            parent[row["name"]] = row["parent"] or None
            rank[row["name"]] = row["rank"]
    if parent:
        RETRO_PARENT, RANK_OF = parent, rank
    return RETRO_PARENT


def load_erv_class(path: str | Path) -> dict[str, str]:
    """
    Replace the ERV-class map from a data file (``erv_class.tsv``: ``genus, erv_class``).

    Curated + overridable; the literal ``ERV_CLASS`` above is only a fallback default.
    """
    mapping: dict[str, str] = {}
    with Path(path).open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            mapping[row["genus"]] = row["erv_class"]
    if mapping:
        ERV_CLASS.clear()
        ERV_CLASS.update(mapping)
    return ERV_CLASS


# =============================================================================
# 2. Pure taxonomy primitives
# =============================================================================
def ancestors(node: str) -> list[str]:
    """Return ``node`` and all its ancestors up to the root (node first)."""
    chain: list[str] = []
    current: str | None = node
    while current is not None:
        chain.append(current)
        current = RETRO_PARENT.get(current)
    return chain


def lca(genera: set[str]) -> str:
    """
    Lowest common ancestor of a set of taxonomy nodes.

    Unknown labels (not in the taxonomy) are ignored; if nothing is known the
    result is ``UNCLASSIFIED``. A single known node returns itself.
    """
    known = [g for g in genera if g in RETRO_PARENT]
    if not known:
        return UNCLASSIFIED
    # Intersect ancestor chains; the deepest shared node is the LCA.
    common: set[str] | None = None
    for g in known:
        chain = set(ancestors(g))
        common = chain if common is None else (common & chain)
    assert common is not None
    # Deepest = the one whose own ancestor list is longest (furthest from root).
    return max(common, key=lambda n: len(ancestors(n)))


def rank_of(node: str) -> str:
    """Taxonomic rank of a resolved node (``genus``/``subfamily``/``family``)."""
    return RANK_OF.get(node, "none")


def weighted_lca(
    hits: list[tuple[str, float]], top_percent: float = 0.10
) -> tuple[str, float]:
    """
    Bitscore-weighted LCA over a locus's hits (MEGAN top-percent paradigm).

    ``hits`` is ``[(genus, bitscore), ...]``. Only hits whose bitscore is within
    ``top_percent`` of the locus's best bitscore count toward the LCA — so a genus
    supported only by weak hits does **not** drag the assignment up to a higher
    rank. This is the fix for unweighted LCA's over-backoff.

    Returns ``(node, confidence)`` where confidence is the bitscore-mass fraction
    carried by the dominant genus among the retained (top) hits: ~1.0 for a clean
    single-genus call, ~1/k when k genera tie near the top.
    """
    known = [(g, s) for g, s in hits if g in RETRO_PARENT and s > 0]
    if not known:
        return UNCLASSIFIED, 0.0
    best = max(s for _, s in known)
    cutoff = best * (1.0 - top_percent)
    kept = [(g, s) for g, s in known if s >= cutoff]

    mass: dict[str, float] = {}
    for g, s in kept:
        mass[g] = mass.get(g, 0.0) + s
    total = sum(mass.values())
    dominant_fraction = max(mass.values()) / total if total else 0.0

    node = lca(set(mass))
    return node, dominant_fraction


# =============================================================================
# 3. GFF3 evidence parsing
# =============================================================================
_LABEL_RE = re.compile(r"label=([^;\t]+)")
_PROBE_RE = re.compile(r"probe=([^;\t]+)")


def parse_genus_set(label_field: str) -> set[str]:
    """
    Parse a GFF3 ``label=`` value into a set of genus names.

    Values are GFF3-escaped (``%3b`` = ``;``) and joined by ``"; "`` when the
    run used list/concatenate aggregation.
    """
    decoded = label_field.replace("%3b", ";").replace("%3B", ";")
    return {part.strip() for part in decoded.split(";") if part.strip()}


def classify_gff3(gff3_path: str | Path) -> list[dict[str, str]]:
    """
    Resolve every feature of a valid-tier GFF3 to an LCA taxon.

    Returns one record per feature with the input genus set, the resolved node,
    its rank, and (when the node is a genus) its ERV class.
    """
    records: list[dict[str, str]] = []
    with Path(gff3_path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = fields[8]
            label_match = _LABEL_RE.search(attrs)
            probe_match = _PROBE_RE.search(attrs)
            if not label_match:
                continue
            genus_set = parse_genus_set(label_match.group(1))
            node = lca(genus_set)
            records.append(
                {
                    "seqname": fields[0],
                    "start": fields[3],
                    "end": fields[4],
                    "probe": probe_match.group(1) if probe_match else "",
                    "n_genera": str(len(genus_set)),
                    "genus_set": ";".join(sorted(genus_set)),
                    "lca_node": node,
                    "lca_rank": rank_of(node),
                    "erv_class": ERV_CLASS.get(node, ""),
                }
            )
    return records


# =============================================================================
# 4. Summary + CLI
# =============================================================================
def summarise(records: list[dict[str, str]]) -> str:
    """Human-readable summary of resolved ranks and confident-genus calls."""
    total = len(records)
    by_rank = Counter(r["lca_rank"] for r in records)
    confident = Counter(r["lca_node"] for r in records if r["lca_rank"] == "genus")
    lines = [f"loci classified: {total}", "", "resolved rank distribution:"]
    for rank in ("genus", "subfamily", "family", "none"):
        if by_rank.get(rank):
            pct = 100.0 * by_rank[rank] / total if total else 0.0
            lines.append(f"  {rank:10s} {by_rank[rank]:6d}  ({pct:4.1f}%)")
    lines += ["", "confident genus calls:"]
    for genus, n in confident.most_common():
        lines.append(f"  {genus:18s} {n:6d}  [{ERV_CLASS.get(genus, '')}]")
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: taxonomy_lca.py <valid.gff3> [out.csv]", file=sys.stderr)
        return 2
    records = classify_gff3(argv[1])
    print(summarise(records))
    if len(argv) >= 3 and records:
        with Path(argv[2]).open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(records[0].keys()))
            writer.writeheader()
            writer.writerows(records)
        print(f"\nwrote per-locus calls -> {argv[2]}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
