"""
Build the taxonomy hierarchy from NCBI Taxonomy (no hard-coding)
===============================================================

Derives the parent->child hierarchy used by `taxonomy_lca` directly from NCBI
Taxonomy, for whatever genera appear in the reference table. Replaces the
hard-coded parent map: the hierarchy now self-updates from the authoritative
source and covers any taxa the reference contains (probe-agnostic at the taxonomy
level), fixing the silent-staleness failure (e.g. the spuma genus rename).

Output: `<reference_dir>/taxonomy.tsv` with columns ``name, parent, rank``
(root's parent is empty). Loaded by `taxonomy_lca.load_taxonomy`.

ERV class (I/II/III) is NOT emitted here: it is a curated biological grouping
(Jern/Blomberg), not an NCBI rank, so it stays a small explicit map in
`taxonomy_lca`. Only the genuine taxonomic hierarchy is derived.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

from Bio import Entrez

ROOT = "Retroviridae"  # keep the subtree from this node down


def genera_from_reference(ref_csv: Path) -> list[str]:
    with ref_csv.open(encoding="utf-8") as fh:
        return sorted({r["genus"] for r in csv.DictReader(fh)})


def lineage_for(genus: str, email: str) -> list[tuple[str, str]] | None:
    """Return [(name, rank), ...] from root down to ``genus`` per NCBI Taxonomy."""
    Entrez.email = email  # type: ignore[assignment]
    with Entrez.esearch(db="taxonomy", term=f"{genus}[Scientific Name]") as h:  # type: ignore[no-untyped-call]
        ids = Entrez.read(h)["IdList"]  # type: ignore[no-untyped-call]
    if not ids:
        return None
    with Entrez.efetch(db="taxonomy", id=ids[0]) as h:  # type: ignore[no-untyped-call]
        rec = Entrez.read(h)[0]  # type: ignore[no-untyped-call]
    chain = [(n["ScientificName"], n.get("Rank", "no rank")) for n in rec["LineageEx"]]
    chain.append((rec["ScientificName"], rec.get("Rank", "genus")))
    return chain


def build(
    genera: list[str], email: str
) -> tuple[dict[str, str | None], dict[str, str]]:
    parent: dict[str, str | None] = {}
    rank: dict[str, str] = {}
    for genus in genera:
        chain = lineage_for(genus, email)
        if not chain:
            print(f"WARN: no NCBI taxid for {genus}", file=sys.stderr)
            continue
        names = [c[0] for c in chain]
        if ROOT in names:  # trim to the retroviral subtree
            chain = chain[names.index(ROOT) :]
        prev: str | None = None
        for name, rk in chain:
            parent.setdefault(name, prev)
            rank[name] = rk
            prev = name
        time.sleep(0.34)  # polite Entrez rate limit
    return parent, rank


def main(argv: list[str]) -> int:
    ref_csv = (
        Path(argv[1])
        if len(argv) > 1
        else Path("data/taxonomy_dev/reference/retro_reference.csv")
    )
    email = argv[2] if len(argv) > 2 else "retroseek@example.org"
    out = ref_csv.parent / "taxonomy.tsv"

    genera = genera_from_reference(ref_csv)
    parent, rank = build(genera, email)
    if not parent:
        print("ERROR: no taxonomy resolved", file=sys.stderr)
        return 1
    with out.open("w", encoding="utf-8") as fh:
        fh.write("name\tparent\trank\n")
        for name, par in parent.items():
            fh.write(f"{name}\t{par or ''}\t{rank.get(name, 'no rank')}\n")
    print(f"wrote {len(parent)} taxonomy nodes -> {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
