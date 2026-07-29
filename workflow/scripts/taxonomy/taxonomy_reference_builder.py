"""
Taxonomy reference builder
==========================

Builds an *independent*, taxon-balanced retroviral protein reference for
weighted-LCA / placement classification of ERV loci (see ``taxonomy_lca.py`` and
``docs/taxonomy_classification/``).

Independence matters: the classifier must reclassify loci against references that
were NOT the probes used to find them, otherwise the assignment is circular. So we
pull RefSeq proteins per **axis taxon** straight from NCBI, balanced by capping the
number kept per (taxon, gene) - this prevents over-sequenced lineages (e.g. MLV
gammaretroviruses) from dominating the reference and biasing the LCA.

The **axis** (the taxa the reference is built at) is rank-agnostic and declared, not
hard-coded (ADR-008): ``classification.reference_taxa`` (``--taxa``) if set, else the
distinct probeset ``Label`` values (``--probe-csv``), else the retroviral-genus default.
``"{taxon}"[Organism]`` matches NCBI at any rank, so a family (``Bornaviridae``) is as
valid an axis member as a genus (``Lentivirus``).

Output (pinned under ``data/taxonomy_reference/`` by the ``taxonomy_reference`` rule):
* ``retro_reference.faa`` - protein FASTA, header = accession.
* ``retro_reference.csv`` - accession, taxon, gene, defline (the taxonomy map the
  weighted-LCA / blastx step joins hits against).
* ``manifest.yaml`` - provenance: axis taxa queried, per-taxon kept counts, total,
  the content hash of the reference CSV, and the Biopython version that fetched it.

Reuses the project's Biopython/Entrez convention (see ``seq_utils.gb_fetcher``).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import logging
import re
import time
from collections import Counter
from pathlib import Path

import Bio
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

from colored_logging import colored_logging

logger = logging.getLogger(__name__)

# Fallback axis when neither --taxa nor --probe-csv is given (CLI/standalone use):
# the retroviral genera. ICTV renamed the spuma genera, so "Spumaretrovirus"[Organism]
# is empty - use the modern spuma genera for Class III.
DEFAULT_TAXA: tuple[str, ...] = (
    "Alpharetrovirus",
    "Betaretrovirus",
    "Gammaretrovirus",
    "Deltaretrovirus",
    "Epsilonretrovirus",
    "Lentivirus",
    "Bovispumavirus",
    "Felispumavirus",
    "Simiispumavirus",
    "Prosimiispumavirus",
    "Equispumavirus",
)

# Gene assignment from the protein defline (case-insensitive, first match wins).
# Comprehensive: anything unmatched is kept as 'OTHER' (accessory/hypothetical), so
# non-canonical genes like REX/TAX classify by taxon too (probe-agnostic at gene level).
GENE_PATTERNS: tuple[tuple[str, str], ...] = (
    ("POL", r"\b(pol|reverse transcriptase|gag-pol|integrase|RNase ?H)\b"),
    ("ENV", r"\b(env|envelope|surface (glyco)?protein|transmembrane)\b"),
    ("GAG", r"\b(gag|capsid|matrix|nucleocapsid)\b"),
    ("PRO", r"\b(pro|protease|proteinase|dUTPase)\b"),
    ("REX", r"\brex\b"),
    ("TAX", r"\btax\b"),
)
OTHER_GENE = "OTHER"

CAP_PER_TAXON_GENE = 25  # balance: keep at most this many per (taxon, gene)
ESEARCH_RETMAX = 600
EFETCH_BATCH = 100


def resolve_axis(taxa: list[str], probe_csv: Path | None) -> list[str]:
    """Resolve the classification axis (ADR-008, hybrid, rank-agnostic).

    Priority: explicit ``taxa`` (``classification.reference_taxa``) > distinct probeset
    ``Label`` values (``probe_csv``) > the retroviral-genus default. Order-preserving,
    de-duplicated.
    """
    if taxa:
        source = taxa
    elif probe_csv is not None and probe_csv.exists():
        labels = pd.read_csv(probe_csv)["Label"].dropna().astype(str)
        source = labels.tolist()
        logger.info("axis derived from %d probeset Label values", len(source))
    else:
        source = list(DEFAULT_TAXA)
        logger.info("axis: no --taxa/--probe-csv given, using retroviral-genus default")
    seen: dict[str, None] = {}
    for t in (s.strip() for s in source):
        if t:
            seen.setdefault(t, None)
    return list(seen)


def classify_gene(defline: str) -> str:
    """Map a protein defline to a coarse gene; unmatched -> OTHER (kept, not dropped)."""
    low = defline.lower()
    for gene, pattern in GENE_PATTERNS:
        if re.search(pattern, low):
            return gene
    return OTHER_GENE


def fetch_taxon(taxon: str, email: str) -> list[tuple[str, str, str, str]]:
    """
    Fetch RefSeq protein records for one axis taxon (any rank).

    Returns list of (accession, taxon, gene, defline), capped per gene for balance.
    """
    Entrez.email = email  # type: ignore[assignment]
    term = f'"{taxon}"[Organism] AND srcdb_refseq[PROP]'
    with Entrez.esearch(db="protein", term=term, retmax=ESEARCH_RETMAX) as h:  # type: ignore[no-untyped-call]
        ids = Entrez.read(h)["IdList"]  # type: ignore[no-untyped-call]
    if not ids:
        return []

    kept: list[tuple[str, str, str, str]] = []
    per_gene: dict[str, int] = {}
    seqs: dict[str, SeqRecord] = {}
    for start in range(0, len(ids), EFETCH_BATCH):
        batch = ids[start : start + EFETCH_BATCH]
        with Entrez.efetch(  # type: ignore[no-untyped-call]
            db="protein", id=",".join(batch), rettype="fasta", retmode="text"
        ) as h:
            for rec in SeqIO.parse(h, "fasta"):  # type: ignore[no-untyped-call]
                gene = classify_gene(rec.description)  # always set (OTHER catch-all)
                if per_gene.get(gene, 0) >= CAP_PER_TAXON_GENE:
                    continue
                acc = rec.id.split()[0]
                per_gene[gene] = per_gene.get(gene, 0) + 1
                kept.append((acc, taxon, gene, rec.description))
                seqs[acc] = rec
        time.sleep(0.4)  # polite Entrez rate limit

    # stash the sequences on the function for the caller to write
    fetch_taxon.last_seqs = seqs  # type: ignore[attr-defined]
    return kept


def write_manifest(
    rows: list[tuple[str, str, str, str]],
    axis: list[str],
    csv_path: Path,
    manifest: Path,
) -> None:
    """Provenance: axis taxa, per-taxon kept counts, total, CSV hash, Biopython version.

    Lets a result be traced to the exact reference snapshot it was classified
    against (the Entrez fetch is the one non-deterministic build step).
    """
    by_taxon = Counter(taxon for _, taxon, _, _ in rows)
    content_hash = hashlib.md5(csv_path.read_bytes()).hexdigest()
    lines = [
        "# RetroSeek taxonomic-classification reference - provenance manifest",
        f"biopython_version: {Bio.__version__}",
        "query_filter: srcdb_refseq[PROP]",
        f"cap_per_taxon_gene: {CAP_PER_TAXON_GENE}",
        f"total_proteins: {len(rows)}",
        f"reference_csv_md5: {content_hash}",
        "axis_taxa:",
    ]
    lines += [f"  {taxon}: {by_taxon.get(taxon, 0)}" for taxon in axis]
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Build the taxon-comprehensive ERV reference"
    )
    p.add_argument("--email", required=True, help="NCBI Entrez contact email")
    p.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="reference package directory (faa/csv/manifest written here)",
    )
    p.add_argument(
        "--taxa",
        default="",
        help="comma-separated axis taxa (classification.reference_taxa); any rank. "
        "Overrides --probe-csv when non-empty.",
    )
    p.add_argument(
        "--probe-csv",
        type=Path,
        default=None,
        help="probeset CSV; axis defaults to its distinct Label values when --taxa is empty",
    )
    a = p.parse_args(argv)
    colored_logging(log_file_name="taxonomy_reference_builder.txt")
    a.out_dir.mkdir(parents=True, exist_ok=True)

    taxa = [t.strip() for t in a.taxa.split(",") if t.strip()]
    axis = resolve_axis(taxa, a.probe_csv)
    logger.info("axis (%d taxa): %s", len(axis), ", ".join(axis))

    all_rows: list[tuple[str, str, str, str]] = []
    all_seqs: dict[str, SeqRecord] = {}
    for taxon in axis:
        rows = fetch_taxon(taxon, a.email)
        all_rows.extend(rows)
        all_seqs.update(fetch_taxon.last_seqs)  # type: ignore[attr-defined]
        by_gene: dict[str, int] = {}
        for _, _, gene, _ in rows:
            by_gene[gene] = by_gene.get(gene, 0) + 1
        logger.info("%-20s kept %3d  %s", taxon, len(rows), by_gene)

    faa = a.out_dir / "retro_reference.faa"
    meta = a.out_dir / "retro_reference.csv"
    with faa.open("w", encoding="utf-8") as fh:
        for acc, *_ in all_rows:
            rec = all_seqs[acc]
            fh.write(f">{acc}\n{rec.seq!s}\n")
    with meta.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["accession", "taxon", "gene", "defline"])
        writer.writerows(all_rows)
    write_manifest(all_rows, axis, meta, a.out_dir / "manifest.yaml")
    logger.info("wrote %d reference proteins -> %s", len(all_rows), faa)
    logger.info("wrote taxonomy map -> %s", meta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
