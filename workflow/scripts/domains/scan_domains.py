# =============================================================================
# scan_domains.py
# =============================================================================
# Protein-domain evidence for every catalogued locus, in BOTH tiers, by one
# identical procedure.
#
# The problem it solves: domains previously reached the catalog only through
# `gt ltrdigest -hmms`, which searches inside LTRharvest elements. Orphan loci sit
# outside every element, so they were hard-defaulted to `non_domain` - 99.04% of
# `non_domain` rows meant "never assessed", indistinguishable from "assessed and
# empty". This scans both tiers and records which it was.
#
# Method symmetry is the whole point, so every axis is held equal:
#   region      the catalogued locus span, on both sides. An LTR element averages
#               6,919 bp against an orphan locus's 681 bp, and a larger search
#               space wins more hits by construction; the locus span is the only
#               region both tiers possess.
#   library     one curated Pfam subset (see subset_pfam.py).
#   threshold   `--cut_ga`, Pfam's per-family curated bit-score cutoffs. Bit
#               scores do not depend on how the search was framed, so these are
#               comparable in a way E-values are not.
#   frames      all six.
#
# Results are keyed on the SET of distinct families per locus, never on hit
# counts. `gt ltrdigest` chains fragments of one model into a single feature and
# `hmmsearch` does not, so any count-based statistic would differ between the two
# for reasons that have nothing to do with biology.
#
# Locus identity is `{seqname}|{parent}`, the natural key shared with the
# classifier. The catalog's own `L{i}` ids are positional and must not be used to
# join across processes.
# =============================================================================

from __future__ import annotations

import argparse
import logging
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd
from Bio.Seq import Seq

# Locus grouping MUST match the classifier exactly or the join silently drops
# rows, so the derivation is imported rather than reimplemented. The pure class
# semantics live in domain_classes, which the classifier also imports; keeping
# them there is what stops this import from becoming a cycle.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "taxonomy"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from taxonomy_classify_loci import build_loci, parse_valid_full

from colored_logging import colored_logging

logger = logging.getLogger(__name__)

HMMSEARCH = "hmmsearch"
_EXTRACT_R = (
    Path(__file__).resolve().parent.parent / "taxonomy" / "extract_region_fasta.R"
)
SIX_FRAMES = (1, 2, 3, -1, -2, -3)


# ----------------------------------------------------------- translation
def six_frame_translations(dna: str) -> dict[int, str]:
    """Translate `dna` in all six frames, keyed by frame (1..3, -1..-3).

    Stop codons become `X` rather than staying `*`, matching
    taxonomy_classify_loci.py. This is not cosmetic: measured on the Homo
    LTR-flanked set, keeping `*` covered 1,946 loci against 2,192 with `X`, and
    1,557 retroviral-diagnostic loci against 1,866. `*` is a non-residue that an
    alignment cannot cross, so a degraded ERV's domain is broken at every stop;
    `X` is the neutral unknown, letting the domain score as one unit.

    Splitting into ORFs at stops instead (the esl-translate convention) was also
    benchmarked and is worse for this data: it recovers a wider family repertoire
    but loses the most degraded loci, because a 102 aa domain broken by one stop
    becomes two ~50 aa fragments that each fail the gathering threshold. See
    ADR-016 and the Phase 1 bench.

    Sequences shorter than a codon yield empty strings rather than raising, which
    keeps a degenerate locus from aborting a whole genome's scan.
    """
    out: dict[int, str] = {}
    for frame in SIX_FRAMES:
        s = Seq(dna)
        if frame < 0:
            s = s.reverse_complement()  # type: ignore[no-untyped-call]
        s = s[abs(frame) - 1 :]
        s = s[: len(s) - (len(s) % 3)]
        out[frame] = str(s.translate()).replace("*", "X")  # type: ignore[no-untyped-call]
    return out


def write_query_fasta(sequences: dict[str, str], path: Path) -> int:
    """Write six-frame translations of `sequences` as `{locus}|f{frame}` records.

    Empty translations are skipped: HMMER rejects zero-length records.
    """
    written = 0
    with path.open("w", encoding="utf-8") as fh:
        for locus_id, dna in sequences.items():
            for frame, protein in six_frame_translations(dna).items():
                if not protein:
                    continue
                fh.write(f">{locus_id}|f{frame}\n{protein}\n")
                written += 1
    return written


def locus_of(target: str) -> str:
    """Strip the `|f{frame}` suffix. Split from the right: locus ids contain `|`."""
    return target.rsplit("|", 1)[0]


# --------------------------------------------------------------- parsing
def parse_domtblout(path: Path) -> list[dict[str, Any]]:
    """Read `hmmsearch --domtblout` rows.

    Field 1 is the TARGET sequence and field 4 the QUERY model; field 5 is the
    model accession. Reversing target and query is the classic mistake with this
    format and produces plausible-looking nonsense rather than an error.
    """
    hits: list[dict[str, Any]] = []
    with path.open(encoding="utf-8") as fh:
        for raw in fh:
            if raw.startswith("#"):
                continue
            f = raw.split()
            if len(f) < 14:
                continue
            target = f[0]
            hits.append(
                {
                    "locus_id": locus_of(target),
                    "frame": int(target.rsplit("|f", 1)[1]),
                    "pfam_name": f[3],
                    "pfam_acc": f[4].split(".")[0],
                    "evalue": float(f[12]),  # i-Evalue, independent of other doms
                    "bitscore": float(f[13]),  # per-domain bit score
                }
            )
    return hits


# --------------------------------------------------------- orchestration
def run(cmd: list[str]) -> None:
    logger.info("running: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def write_locus_bed(loci: list[dict[str, Any]], bed: Path) -> None:
    """One BED row per locus, spanning the whole catalogued locus.

    Distinct from the classifier's `write_region_bed`, which emits one row per
    (locus, gene): the scan deliberately measures the locus, not its gene parts.
    Strand is forced to `+` because all six frames are translated anyway.
    """
    with bed.open("w", encoding="utf-8") as fh:
        for lc in loci:
            key = f"{lc['seqname']}|{lc['parent']}"
            fh.write(
                f"{lc['seqname']}\t{max(0, int(lc['start']) - 1)}\t{int(lc['end'])}"
                f"\t{key}\t.\t+\n"
            )


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def scan_regions(
    bed: Path,
    genome: Path,
    hmm: Path,
    workdir: Path,
    source: str,
    threads: int,
) -> list[dict[str, Any]]:
    """Extract, translate and search the regions in `bed`. Returns raw hits.

    Shared by the locus scan and the element scan. Keeping one body is what makes
    the two comparable: they differ only in which regions the BED names.
    """
    fna = workdir / "regions.fna"
    run(
        [
            "Rscript",
            str(_EXTRACT_R),
            "--genome",
            str(genome),
            "--bed",
            str(bed),
            "--out",
            str(fna),
        ]
    )

    queries = workdir / "queries.faa"
    n_queries = write_query_fasta(read_fasta(fna), queries)
    logger.info("%s: %d translated query frames", source, n_queries)

    domtbl = workdir / "hits.domtbl"
    if n_queries:
        run(
            [
                HMMSEARCH,
                "--cut_ga",
                "--cpu",
                str(threads),
                "-o",
                "/dev/null",
                "--domtblout",
                str(domtbl),
                str(hmm),
                str(queries),
            ]
        )
    else:
        domtbl.write_text("# no queries\n")

    hits = parse_domtblout(domtbl)
    logger.info("%s: %d domain hits at GA", source, len(hits))
    for hit in hits:
        hit["source"] = source
    return hits


_HIT_COLUMNS = [
    "locus_id",
    "source",
    "pfam_acc",
    "pfam_name",
    "bitscore",
    "evalue",
    "frame",
]


def scan(
    gff3: Path,
    genome: Path,
    hmm: Path,
    workdir: Path,
    source: str,
    threads: int,
) -> pd.DataFrame:
    """Scan every catalogued locus in `gff3`, one row per (locus, domain hit)."""
    workdir.mkdir(parents=True, exist_ok=True)
    loci = build_loci(parse_valid_full(gff3))
    logger.info("%s: %d loci to scan", source, len(loci))
    bed = workdir / "regions.bed"
    write_locus_bed(loci, bed)
    hits = scan_regions(bed, genome, hmm, workdir, source, threads)
    frame = pd.DataFrame(hits, columns=_HIT_COLUMNS)
    # Loci that produced no hit still need a row in the provenance record, so the
    # classifier can tell `non_domain` (scanned, empty) from `not_scanned`.
    frame.attrs["scanned"] = {f"{lc['seqname']}|{lc['parent']}" for lc in loci}
    return frame


def _run_loci(args: argparse.Namespace) -> None:
    """Locus mode: the symmetric cross-tier scan that feeds the catalog columns."""
    frames, scanned = [], set()
    for gff3, source in (
        (args.valid_gff3, "ltr-flanked"),
        (args.orphan_gff3, "orphan"),
    ):
        df = scan(
            gff3, args.genome, args.hmm, args.workdir / source, source, args.threads
        )
        scanned |= df.attrs["scanned"]
        frames.append(df)

    out = pd.concat(frames, ignore_index=True)
    args.out_parquet.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(args.out_parquet, index=False)
    out.to_csv(args.out_csv, index=False)
    # The scanned-locus roster is what makes "no domains" distinguishable from
    # "never looked"; without it the old ambiguity returns.
    args.out_scanned.write_text("\n".join(sorted(scanned)) + "\n", encoding="utf-8")
    logger.info("wrote %d domain hits over %d scanned loci", len(out), len(scanned))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Scan catalogued loci for Pfam domains"
    )
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--hmm", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--out-parquet", type=Path, required=True)
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--valid-gff3", type=Path, required=True)
    parser.add_argument("--orphan-gff3", type=Path, required=True)
    parser.add_argument("--out-scanned", type=Path, required=True)
    args = parser.parse_args()
    colored_logging(log_file_name=f"domain_scan_{args.genome.stem}.txt")

    _run_loci(args)


if __name__ == "__main__":
    main()
