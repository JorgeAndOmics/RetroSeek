"""Build a per-gene reference tree package for phylogenetic placement.

For one gene (e.g. POL, GAG), builds the artifacts EPA-ng/gappa need, from the
gene's subset of the taxon-comprehensive reference:

    gene refs (faa)  --MAFFT-->  <gene>.afa      reference MSA
                     --IQ-TREE-> <gene>.treefile + <gene>.model   tree + substitution model
                     --hmmbuild->  <gene>.hmm    profile (published, NOT used by placement)
    + <gene>.taxon.tsv   tip(accession) -> lineage   (for gappa examine assign)

``<gene>.hmm`` is not on the placement path. Queries are aligned with
``mafft --add --keeplength`` against ``<gene>.afa`` (see
`taxonomy_placement._align_queries`), and EPA-ng takes ``--tree``/``--ref-msa``/
``--query``/``--model`` - it accepts no profile. Both routes enforce the same
invariant EPA-ng requires, that a query occupies exactly the reference's
columns; ``hmmalign`` against this profile was the design's other option and
MAFFT was chosen instead. The profile is still built and published because it
costs under a second on a ~64-sequence alignment and lets you hmmsearch your own
sequences against the reference without rebuilding anything. Do not infer from
its presence that HMMER participates in placement.

Also reports alignment quality (n seqs, columns, % gaps, mean pairwise identity) so
poorly-aligning markers (e.g. ENV) are flagged as low-confidence trees.

Tools (mafft, iqtree, raxml-ng, hmmbuild) resolve from PATH in the `RetroSeek` env.
Pinned under <ref_dir>/trees/. See docs/taxonomy_classification/.
"""

from __future__ import annotations

import argparse
import csv
import logging
from itertools import combinations
from pathlib import Path

import taxonomy_lca as tlca
from Bio import SeqIO

from external import run_tool
from log import OK, PipelineError, job_logging, run_main

logger = logging.getLogger(__name__)

MAFFT = "mafft"  # all tools resolved from PATH (the RetroSeek conda env)
IQTREE = "iqtree"
RAXML = "raxml-ng"
HMMBUILD = "hmmbuild"


def gene_subset(
    ref_csv: Path, ref_faa: Path, gene: str, out_faa: Path
) -> dict[str, str]:
    """Write the gene's reference proteins to out_faa; return accession->taxon."""
    wanted = {}
    with ref_csv.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["gene"] == gene:
                wanted[r["accession"]] = r["taxon"]
    n = 0
    with out_faa.open("w", encoding="utf-8") as out:
        for rec in SeqIO.parse(str(ref_faa), "fasta"):  # type: ignore[no-untyped-call]
            acc = rec.id.split()[0]
            if acc in wanted:
                out.write(f">{acc}\n{rec.seq!s}\n")
                n += 1
    if n < 4:
        raise PipelineError(
            f"only {n} {gene} reference proteins, too few to build a tree",
            hint="widen classification.reference_taxa or drop the gene from placement_genes",
        )
    return wanted


def _pair_identity(a: str, b: str) -> float | None:
    """Fraction identical over the columns where neither sequence has a gap.

    None when the two share no such column.
    """
    cols = [(x, y) for x, y in zip(a, b, strict=False) if x != "-" and y != "-"]
    if not cols:
        return None
    return sum(x == y for x, y in cols) / len(cols)


def _mean_pairwise_identity(seqs: list[str], max_pairs: int = 40) -> float:
    """Mean % identity over the first ``max_pairs`` pairs that share a column.

    Pairs are taken in alignment order, (0,1), (0,2), ..., (1,2), ..., so the
    value depends on the ORDER of the sequences, not only on the alignment: the
    same sequences in another order sample other pairs. Kept that way so the
    logged number stays comparable across runs. Re-measuring it on a reordered
    subset and calling the documented value stale is a known trap.
    """
    # A plain running total on purpose: from Python 3.12 sum() compensates float
    # error and would move the last bits of the logged value.
    total, pairs = 0.0, 0
    for a, b in combinations(seqs, 2):
        identity = _pair_identity(a, b)
        if identity is None:
            continue
        total += identity
        pairs += 1
        if pairs == max_pairs:
            break
    return 100.0 * total / pairs if pairs else 0.0


def alignment_quality(afa: Path) -> str:
    """N seqs, columns, % gap, mean pairwise %identity (sampled) - a quality flag."""
    seqs = [str(r.seq) for r in SeqIO.parse(str(afa), "fasta")]  # type: ignore[no-untyped-call]
    ncol = len(seqs[0]) if seqs else 0
    gap = sum(s.count("-") for s in seqs) / max(1, len(seqs) * ncol)
    mpi = _mean_pairwise_identity(seqs)
    flag = "OK" if mpi >= 25 else "LOW (tree may be unreliable)"
    return f"n={len(seqs)} cols={ncol} gap={100 * gap:.0f}% mean_pident={mpi:.1f}% -> {flag}"


def taxon_map(acc_taxon: dict[str, str], out_tsv: Path) -> None:
    """tip(accession) -> 'root;...;taxon' lineage for gappa (any rank; ADR-008)."""
    with out_tsv.open("w", encoding="utf-8") as fh:
        for acc, taxon in acc_taxon.items():
            lineage = ";".join(reversed(tlca.ancestors(taxon)))  # root..taxon
            fh.write(f"{acc}\t{lineage}\n")


def main(argv: list[str] | None = None) -> None:
    """Command-line entry: build the placement tree package for one gene."""
    p = argparse.ArgumentParser(description="Build a per-gene placement tree package")
    p.add_argument("gene", help="marker gene, e.g. POL or GAG")
    p.add_argument(
        "--ref-dir",
        type=Path,
        required=True,
        help="reference package dir (reads retro_reference.{faa,csv}, taxonomy.tsv)",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=0,
        help="RNG seed for IQ-TREE + raxml-ng (reproducible trees); from parameters.seed",
    )
    p.add_argument("--threads", type=int, default=2, help="tree-building threads")
    p.add_argument("--log", type=Path, help="job log (the Snakemake log: path)")
    a = p.parse_args(argv)
    gene = a.gene.upper()
    job_logging(a.log, "taxonomy_reference_trees")
    ref_dir = a.ref_dir
    trees = ref_dir / "trees"
    trees.mkdir(parents=True, exist_ok=True)

    tax = ref_dir / "taxonomy.tsv"
    if tax.exists():
        tlca.load_taxonomy(tax)

    gene_faa = trees / f"{gene}.faa"
    acc_taxon = gene_subset(
        ref_dir / "retro_reference.csv", ref_dir / "retro_reference.faa", gene, gene_faa
    )

    afa = trees / f"{gene}.afa"
    run_tool(
        [MAFFT, "--maxiterate", "1000", "--localpair", "--anysymbol", str(gene_faa)],
        stdout=afa.open("w", encoding="utf-8"),
    )
    logger.info("[%s] alignment: %s", gene, alignment_quality(afa))

    prefix = trees / gene
    # fixed standard protein model (LG+F+G4) - skips slow ModelFinder; gives the topology.
    # -seed makes the search reproducible (UFBoot resampling + NNI tie-breaks).
    run_tool(
        [
            IQTREE,
            "-s",
            str(afa),
            "-m",
            "LG+F+G4",
            "-B",
            "1000",
            "-T",
            str(a.threads),
            "-seed",
            str(a.seed),
            "--prefix",
            str(prefix),
            "-redo",
        ]
    )
    # raxml-ng --evaluate: optimise model params + branch lengths on the fixed topology ->
    # <gene>.raxml.bestModel + .bestTree, which EPA-ng uses for calibrated placement.
    # raxml-ng appends ".raxml.<suffix>" to --prefix, so prefix=<gene> -> <gene>.raxml.bestModel
    # --redo overwrites stale <gene>.raxml.* from a prior build; without it raxml-ng aborts
    # on any rebuild ("file already exists"), unlike the idempotent iqtree (-redo) call.
    run_tool(
        [
            RAXML,
            "--evaluate",
            "--redo",
            "--msa",
            str(afa),
            "--tree",
            str(prefix) + ".treefile",
            "--model",
            "LG+F+G4",
            "--prefix",
            str(prefix),
            "--seed",
            str(a.seed),
            "--force",
            "perf_threads,msa",
            "--threads",
            str(a.threads),
        ]
    )
    run_tool([HMMBUILD, "--amino", str(trees / f"{gene}.hmm"), str(afa)])
    taxon_map(acc_taxon, trees / f"{gene}.taxon.tsv")
    logger.log(OK, "%s tree package written to %s", gene, trees)


if __name__ == "__main__":
    run_main(main)
