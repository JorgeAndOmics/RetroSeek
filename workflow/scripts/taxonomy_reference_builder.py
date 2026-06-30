"""
Taxonomy reference builder
==========================

Builds an *independent*, genus-balanced retroviral protein reference for
weighted-LCA classification of ERV loci (see ``taxonomy_lca.py`` and
``docs/taxonomy_classification/``).

Independence matters: the classifier must reclassify loci against references that
were NOT the probes used to find them, otherwise the assignment is circular. So we
pull RefSeq proteins per retroviral genus straight from NCBI, balanced by capping
the number kept per (genus, gene) — this prevents over-sequenced lineages (e.g.
MLV gammaretroviruses) from dominating the reference and biasing the LCA.

Output (pinned under ``data/taxonomy_reference/`` by the ``taxonomy_reference`` rule):
* ``retro_reference.faa`` — protein FASTA, header = accession.
* ``retro_reference.csv`` — accession, genus, gene, defline (the taxonomy map the
  weighted-LCA / blastx step joins hits against).
* ``manifest.yaml`` — provenance: genera queried, per-genus kept counts, total,
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
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

from colored_logging import colored_logging

logger = logging.getLogger(__name__)

# The 7 ICTV retroviral genera — the classification axis (probe-agnostic: these are
# reference genera, independent of whatever probes a run declares).
# Retroviral genera = the classification axis (probe-agnostic: reference genera,
# independent of whatever probes a run declares). ICTV renamed the spuma genera, so
# "Spumaretrovirus"[Organism] is empty — use the modern spuma genera for Class III.
GENERA: tuple[str, ...] = (
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
# Genus-comprehensive: anything unmatched is kept as 'OTHER' (accessory/hypothetical),
# so non-canonical genes like REX/TAX classify by genus too (probe-agnostic at gene level).
GENE_PATTERNS: tuple[tuple[str, str], ...] = (
    ("POL", r"\b(pol|reverse transcriptase|gag-pol|integrase|RNase ?H)\b"),
    ("ENV", r"\b(env|envelope|surface (glyco)?protein|transmembrane)\b"),
    ("GAG", r"\b(gag|capsid|matrix|nucleocapsid)\b"),
    ("PRO", r"\b(pro|protease|proteinase|dUTPase)\b"),
    ("REX", r"\brex\b"),
    ("TAX", r"\btax\b"),
)
OTHER_GENE = "OTHER"

CAP_PER_GENUS_GENE = 25  # balance: keep at most this many per (genus, gene)
ESEARCH_RETMAX = 600
EFETCH_BATCH = 100


def classify_gene(defline: str) -> str:
    """Map a protein defline to a coarse gene; unmatched -> OTHER (kept, not dropped)."""
    low = defline.lower()
    for gene, pattern in GENE_PATTERNS:
        if re.search(pattern, low):
            return gene
    return OTHER_GENE


def fetch_genus(genus: str, email: str) -> list[tuple[str, str, str, str]]:
    """
    Fetch RefSeq protein records for one genus.

    Returns list of (accession, genus, gene, defline), capped per gene for balance.
    """
    Entrez.email = email  # type: ignore[assignment]
    term = f'"{genus}"[Organism] AND srcdb_refseq[PROP]'
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
                if per_gene.get(gene, 0) >= CAP_PER_GENUS_GENE:
                    continue
                acc = rec.id.split()[0]
                per_gene[gene] = per_gene.get(gene, 0) + 1
                kept.append((acc, genus, gene, rec.description))
                seqs[acc] = rec
        time.sleep(0.4)  # polite Entrez rate limit

    # stash the sequences on the function for the caller to write
    fetch_genus.last_seqs = seqs  # type: ignore[attr-defined]
    return kept


def write_manifest(
    rows: list[tuple[str, str, str, str]], csv_path: Path, manifest: Path
) -> None:
    """Provenance: genera, per-genus kept counts, total, CSV hash, Biopython version.

    Lets a result be traced to the exact reference snapshot it was classified
    against (the Entrez fetch is the one non-deterministic build step).
    """
    by_genus = Counter(genus for _, genus, _, _ in rows)
    content_hash = hashlib.md5(csv_path.read_bytes()).hexdigest()
    lines = [
        "# RetroSeek taxonomic-classification reference — provenance manifest",
        f"biopython_version: {Bio.__version__}",
        "query_filter: srcdb_refseq[PROP]",
        f"cap_per_genus_gene: {CAP_PER_GENUS_GENE}",
        f"total_proteins: {len(rows)}",
        f"reference_csv_md5: {content_hash}",
        "genera:",
    ]
    lines += [f"  {genus}: {by_genus.get(genus, 0)}" for genus in GENERA]
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Build the genus-comprehensive ERV reference"
    )
    p.add_argument("--email", required=True, help="NCBI Entrez contact email")
    p.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="reference package directory (faa/csv/manifest written here)",
    )
    a = p.parse_args(argv)
    colored_logging(log_file_name="taxonomy_reference_builder.txt")
    a.out_dir.mkdir(parents=True, exist_ok=True)

    all_rows: list[tuple[str, str, str, str]] = []
    all_seqs: dict[str, SeqRecord] = {}
    for genus in GENERA:
        rows = fetch_genus(genus, a.email)
        all_rows.extend(rows)
        all_seqs.update(fetch_genus.last_seqs)  # type: ignore[attr-defined]
        by_gene: dict[str, int] = {}
        for _, _, gene, _ in rows:
            by_gene[gene] = by_gene.get(gene, 0) + 1
        logger.info("%-20s kept %3d  %s", genus, len(rows), by_gene)

    faa = a.out_dir / "retro_reference.faa"
    meta = a.out_dir / "retro_reference.csv"
    with faa.open("w", encoding="utf-8") as fh:
        for acc, *_ in all_rows:
            rec = all_seqs[acc]
            fh.write(f">{acc}\n{rec.seq!s}\n")
    with meta.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["accession", "genus", "gene", "defline"])
        writer.writerows(all_rows)
    write_manifest(all_rows, meta, a.out_dir / "manifest.yaml")
    logger.info("wrote %d reference proteins -> %s", len(all_rows), faa)
    logger.info("wrote taxonomy map -> %s", meta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
