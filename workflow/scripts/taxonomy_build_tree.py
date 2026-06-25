"""
Build a per-gene reference tree package for phylogenetic placement
=================================================================

For one gene (e.g. POL, GAG), builds the artifacts EPA-ng/gappa need, from the
gene's subset of the genus-comprehensive reference:

    gene refs (faa)  --MAFFT-->  <gene>.afa      reference MSA
                     --IQ-TREE-> <gene>.treefile + <gene>.model   tree + substitution model
                     --hmmbuild->  <gene>.hmm    profile (for aligning queries at run time)
    + <gene>.taxon.tsv   tip(accession) -> lineage   (for gappa examine assign)

Also reports alignment quality (n seqs, columns, % gaps, mean pairwise identity) so
poorly-aligning markers (e.g. ENV) are flagged as low-confidence trees.

Tools (mafft, iqtree, raxml-ng, hmmbuild) resolve from PATH in the `RetroSeek` env.
Pinned under <ref_dir>/trees/. See docs/taxonomy_classification/.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path
from typing import Any

from Bio import SeqIO

import taxonomy_lca as tlca

MAFFT = "mafft"  # all tools resolved from PATH (the RetroSeek conda env)
IQTREE = "iqtree"
RAXML = "raxml-ng"
HMMBUILD = "hmmbuild"


def run(cmd: list[str], stdout: Any = None) -> subprocess.CompletedProcess[str]:
    res = subprocess.run(
        cmd,
        stdout=stdout or subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if res.returncode != 0:
        sys.stderr.write((res.stderr or "")[-3000:])
        raise SystemExit(f"command failed: {' '.join(cmd[:3])}…")
    return res


def gene_subset(
    ref_csv: Path, ref_faa: Path, gene: str, out_faa: Path
) -> dict[str, str]:
    """Write the gene's reference proteins to out_faa; return accession->genus."""
    wanted = {}
    with ref_csv.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["gene"] == gene:
                wanted[r["accession"]] = r["genus"]
    n = 0
    with out_faa.open("w", encoding="utf-8") as out:
        for rec in SeqIO.parse(str(ref_faa), "fasta"):  # type: ignore[no-untyped-call]
            acc = rec.id.split()[0]
            if acc in wanted:
                out.write(f">{acc}\n{rec.seq!s}\n")
                n += 1
    if n < 4:
        raise SystemExit(f"only {n} {gene} references — too few to build a tree")
    return wanted


def alignment_quality(afa: Path) -> str:
    """n seqs, columns, % gap, mean pairwise %identity (sampled) — a quality flag."""
    seqs = [str(r.seq) for r in SeqIO.parse(str(afa), "fasta")]  # type: ignore[no-untyped-call]
    ncol = len(seqs[0]) if seqs else 0
    gap = sum(s.count("-") for s in seqs) / max(1, len(seqs) * ncol)
    # mean pairwise identity over non-gap columns (sample up to 40 pairs)
    pairs, idents = 0, 0.0
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            a, b = seqs[i], seqs[j]
            cols = [(x, y) for x, y in zip(a, b, strict=False) if x != "-" and y != "-"]
            if cols:
                idents += sum(x == y for x, y in cols) / len(cols)
                pairs += 1
            if pairs >= 40:
                break
        if pairs >= 40:
            break
    mpi = 100.0 * idents / pairs if pairs else 0.0
    flag = "OK" if mpi >= 25 else "LOW (tree may be unreliable)"
    return f"n={len(seqs)} cols={ncol} gap={100 * gap:.0f}% mean_pident={mpi:.1f}% -> {flag}"


def taxon_map(acc_genus: dict[str, str], out_tsv: Path) -> None:
    """tip(accession) -> 'Family;Subfamily;Genus' lineage for gappa."""
    with out_tsv.open("w", encoding="utf-8") as fh:
        for acc, genus in acc_genus.items():
            lineage = ";".join(reversed(tlca.ancestors(genus)))  # root..genus
            fh.write(f"{acc}\t{lineage}\n")


def main(argv: list[str] | None = None) -> int:
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
    a = p.parse_args(argv)
    gene = a.gene.upper()
    ref_dir = a.ref_dir
    trees = ref_dir / "trees"
    trees.mkdir(parents=True, exist_ok=True)

    tax = ref_dir / "taxonomy.tsv"
    if tax.exists():
        tlca.load_taxonomy(tax)

    gene_faa = trees / f"{gene}.faa"
    acc_genus = gene_subset(
        ref_dir / "retro_reference.csv", ref_dir / "retro_reference.faa", gene, gene_faa
    )

    afa = trees / f"{gene}.afa"
    run(
        [MAFFT, "--maxiterate", "1000", "--localpair", "--anysymbol", str(gene_faa)],
        stdout=afa.open("w", encoding="utf-8"),
    )
    print(f"[{gene}] alignment: {alignment_quality(afa)}", file=sys.stderr)

    prefix = trees / gene
    # fixed standard protein model (LG+F+G4) — skips slow ModelFinder; gives the topology.
    # -seed makes the search reproducible (UFBoot resampling + NNI tie-breaks).
    run(
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
    run(
        [
            RAXML,
            "--evaluate",
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
    run([HMMBUILD, "--amino", str(trees / f"{gene}.hmm"), str(afa)])
    taxon_map(acc_genus, trees / f"{gene}.taxon.tsv")
    print(f"[{gene}] tree package written -> {trees}/{gene}.*", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
