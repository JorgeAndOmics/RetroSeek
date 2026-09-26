"""Build the taxonomy hierarchy from NCBI Taxonomy (no hard-coding).

Derives the parent->child hierarchy used by `taxonomy_lca` directly from NCBI
Taxonomy, for whatever axis taxa appear in the reference table (any rank; ADR-008).
Replaces the hard-coded parent map: the hierarchy now self-updates from the
authoritative source and covers any taxa the reference contains (probe-agnostic at
the taxonomy level), fixing the silent-staleness failure (e.g. the spuma genus rename).

The full NCBI lineage of each axis taxon is kept (no subtree trim), so mixed-rank axes
- a retroviral genus and a non-retroviral family - share one hierarchy and
`taxonomy_lca.lca` can find their true common ancestor.

Output: `<reference_dir>/taxonomy.tsv` with columns ``name, parent, rank``
(root's parent is empty). Loaded by `taxonomy_lca.load_taxonomy`.

ERV class (I/II/III) is NOT emitted here: it is a curated biological grouping
(Jern/Blomberg), not an NCBI rank, so it stays a small explicit map in
`taxonomy_lca`. Only the genuine taxonomic hierarchy is derived.
"""

from __future__ import annotations

import argparse
import csv
import logging
import time
from pathlib import Path

from Bio import Entrez

from log import OK, PipelineError, job_logging, run_main

logger = logging.getLogger(__name__)


def taxa_from_reference(ref_csv: Path) -> list[str]:
    """The distinct axis taxa in the reference CSV's `taxon` column, sorted."""
    with ref_csv.open(encoding="utf-8") as fh:
        return sorted({r["taxon"] for r in csv.DictReader(fh)})


def lineage_for(taxon: str, email: str) -> list[tuple[str, str]] | None:
    """Return [(name, rank), ...] from root down to ``taxon`` per NCBI Taxonomy."""
    Entrez.email = email  # type: ignore[assignment]
    with Entrez.esearch(db="taxonomy", term=f"{taxon}[Scientific Name]") as h:  # type: ignore[no-untyped-call]
        ids = Entrez.read(h)["IdList"]  # type: ignore[no-untyped-call]
    if not ids:
        return None
    with Entrez.efetch(db="taxonomy", id=ids[0]) as h:  # type: ignore[no-untyped-call]
        rec = Entrez.read(h)[0]  # type: ignore[no-untyped-call]
    chain = [(n["ScientificName"], n.get("Rank", "no rank")) for n in rec["LineageEx"]]
    chain.append((rec["ScientificName"], rec.get("Rank", "no rank")))
    return chain


def build(taxa: list[str], email: str) -> tuple[dict[str, str | None], dict[str, str]]:
    """Merge the NCBI lineages of `taxa` into ({name: parent}, {name: rank}).

    The root's parent is None. A taxon NCBI does not know is skipped with a
    warning.
    """
    # Full NCBI lineage of each axis taxon is kept (no subtree trim; ADR-008), so
    # `taxonomy_lca.lca` can resolve the true common ancestor of a mixed-rank axis.
    parent: dict[str, str | None] = {}
    rank: dict[str, str] = {}
    for taxon in taxa:
        chain = lineage_for(taxon, email)
        if not chain:
            logger.warning("no NCBI taxid for %s", taxon)
            continue
        prev: str | None = None
        for name, rk in chain:
            parent.setdefault(name, prev)
            rank[name] = rk
            prev = name
        time.sleep(0.34)  # polite Entrez rate limit
    return parent, rank


def main(argv: list[str] | None = None) -> None:
    """Command-line entry: write taxonomy.tsv beside the reference CSV.

    Stops with a PipelineError when no taxon of the reference resolves in NCBI.
    """
    parser = argparse.ArgumentParser(
        description="Build taxonomy.tsv from NCBI Taxonomy."
    )
    parser.add_argument(
        "ref_csv", type=Path, help="retro_reference.csv of the reference"
    )
    parser.add_argument("email", help="NCBI Entrez contact email")
    parser.add_argument("--log", type=Path, help="job log (the Snakemake log: path)")
    args = parser.parse_args(argv)
    job_logging(args.log, "taxonomy_reference")
    out = args.ref_csv.parent / "taxonomy.tsv"

    taxa = taxa_from_reference(args.ref_csv)
    parent, rank = build(taxa, args.email)
    if not parent:
        raise PipelineError(
            "no taxon of the reference resolved in NCBI Taxonomy",
            hint="check the network and execution.entrez_email, then rerun",
        )
    with out.open("w", encoding="utf-8") as fh:
        fh.write("name\tparent\trank\n")
        for name, par in parent.items():
            fh.write(f"{name}\t{par or ''}\t{rank.get(name, 'no rank')}\n")
    logger.log(OK, "%s taxonomy nodes", f"{len(parent):,}")


if __name__ == "__main__":
    run_main(main)
