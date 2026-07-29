"""
Phylogenetic placement of marker sequences (POL/GAG) onto a reference tree
=========================================================================

`place()` takes a batch of per-locus marker proteins for one gene plus that gene's
pinned tree package, and returns a taxon call (+ rank + confidence + taxopath) per
locus by:

    queries --mafft --add --keeplength--> aligned into the reference's exact columns
            --EPA-ng--> .jplace (which branch each query attaches to, with weights)
            --gappa examine assign--> taxopath (lineage) + confidence per query

This is the placement branch of the classifier dispatcher; non-placement genes use
`taxonomy_lca.weighted_lca`. Tools (mafft, epa-ng, gappa) resolve from PATH in the
`RetroSeek` env. See docs/taxonomy_classification/.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Any

import taxonomy_lca as tlca
from Bio import SeqIO

MAFFT = "mafft"  # all tools resolved from PATH (the RetroSeek conda env)
EPA_NG = "epa-ng"
GAPPA = "gappa"


def _run(cmd: list[str], stdout: Any = None) -> subprocess.CompletedProcess[str]:
    res = subprocess.run(
        cmd,
        stdout=stdout or subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if res.returncode != 0:
        sys.stderr.write((res.stderr or "")[-3000:])
        raise SystemExit(f"placement command failed: {' '.join(cmd[:3])}...")
    return res


def _read_model(iqtree_file: Path, default: str = "LG+F+G4") -> str:
    """Read the substitution model string from an IQ-TREE .iqtree report."""
    if iqtree_file.exists():
        for line in iqtree_file.read_text(encoding="utf-8").splitlines():
            if "Model of substitution:" in line:
                return line.split(":", 1)[1].strip()
    return default


# Standard 20 amino acids. A post-alignment query row must overlap the reference
# columns at enough informative sites to be placeable. epa-ng ABORTS THE WHOLE
# RUN on a query it considers to have "no non-gap sites", and it fires this even
# for a query aligning at a single residue (observed: an all-gap-but-one-'F' POL
# fragment). Gaps (-, .) and fully-ambiguous residues (X, *, ?) are all
# non-informative. Below the threshold a query carries no phylogenetic signal
# anyway, so it is dropped from placement and falls back to weighted-LCA.
_STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
_MIN_PLACEMENT_SITES = 10  # min standard-AA columns for a meaningful placement


def _informative_site_count(seq: str) -> int:
    """Number of standard-amino-acid columns in a (possibly aligned) sequence."""
    return sum(ch in _STANDARD_AA for ch in seq.upper())


def _is_placeable(seq: str) -> bool:
    """True if the aligned query has enough informative sites to place (>= min)."""
    return _informative_site_count(seq) >= _MIN_PLACEMENT_SITES


def _align_queries(
    queries: dict[str, str], ref_afa: Path, workdir: Path
) -> Path | None:
    """Align query markers into the reference's exact columns (mafft --add --keeplength)."""
    q_faa = workdir / "queries.faa"
    with q_faa.open("w", encoding="utf-8") as fh:
        for qid, seq in queries.items():
            fh.write(f">{qid}\n{seq}\n")
    combined = workdir / "combined.afa"
    _run(
        [MAFFT, "--add", str(q_faa), "--keeplength", "--anysymbol", str(ref_afa)],
        stdout=combined.open("w", encoding="utf-8"),
    )
    # keep only the query rows (same columns as ref); DROP rows with no informative
    # residues (all-gap OR all-X) - EPA-ng aborts the whole run on a query with no
    # non-gap sites, so one degenerate fragment must not take down the placement.
    q_aln = workdir / "query_aln.afa"
    qids = set(queries)
    kept = 0
    with q_aln.open("w", encoding="utf-8") as out:
        for rec in SeqIO.parse(str(combined), "fasta"):  # type: ignore[no-untyped-call]
            if rec.id.split()[0] in qids and _is_placeable(str(rec.seq)):
                out.write(f">{rec.id.split()[0]}\n{rec.seq!s}\n")
                kept += 1
    return q_aln if kept else None


def _parse_gappa(per_query_tsv: Path) -> dict[str, tuple[str, float]]:
    """Parse gappa examine assign output -> query_id -> (taxopath, confidence)."""
    out: dict[str, tuple[str, float]] = {}
    if not per_query_tsv.exists():
        return out
    with per_query_tsv.open(encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        name_i = idx.get("name", 0)
        path_i = idx.get("taxopath", len(header) - 1)
        # confidence: prefer aLWR, then LWR, else 0
        conf_i = idx.get("aLWR", idx.get("LWR"))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= path_i:
                continue
            conf = 0.0
            if conf_i is not None and conf_i < len(f):
                try:
                    conf = float(f[conf_i])
                except ValueError:
                    conf = 0.0
            out[f[name_i]] = (f[path_i], conf)
    return out


def place(
    queries: dict[str, str], ref_dir: Path, gene: str, workdir: Path
) -> dict[str, dict[str, str]]:
    """
    Place query markers for `gene`; return query_id -> {taxon_call, rank, confidence}.

    taxon_call = deepest node of the assigned taxopath that is in the taxonomy
    (an axis taxon if resolved, else a higher rank); rank via taxonomy_lca.rank_of.
    The caller accepts the placement only when taxon_call is an axis member (ADR-008).
    """
    workdir.mkdir(parents=True, exist_ok=True)
    trees = ref_dir / "trees"
    ref_afa = trees / f"{gene}.afa"
    treefile = trees / f"{gene}.treefile"
    taxon_tsv = trees / f"{gene}.taxon.tsv"
    if not (ref_afa.exists() and treefile.exists() and taxon_tsv.exists()):
        return {}  # no tree package -> caller falls back to LCA

    # prefer the raxml-ng-optimised model + tree (calibrated likelihoods -> less over-backoff);
    # fall back to the iqtree tree + bare model string if raxml --evaluate wasn't run.
    best_model = trees / f"{gene}.raxml.bestModel"
    best_tree = trees / f"{gene}.raxml.bestTree"
    use_tree = best_tree if best_tree.exists() else treefile
    use_model = (
        str(best_model)
        if best_model.exists()
        else _read_model(trees / f"{gene}.iqtree")
    )
    q_aln = _align_queries(queries, ref_afa, workdir)
    if q_aln is None:  # all queries were all-gap -> nothing to place
        return {}
    _run(
        [
            EPA_NG,
            "--redo",
            "--tree",
            str(use_tree),
            "--ref-msa",
            str(ref_afa),
            "--query",
            str(q_aln),
            "--model",
            use_model,
            "--outdir",
            str(workdir),
        ]
    )
    jplace = workdir / "epa_result.jplace"
    _run(
        [
            GAPPA,
            "examine",
            "assign",
            "--jplace-path",
            str(jplace),
            "--taxon-file",
            str(taxon_tsv),
            "--per-query-results",
            "--best-hit",
            "--out-dir",
            str(workdir),
            "--allow-file-overwriting",
        ]
    )

    assigned = _parse_gappa(workdir / "per_query.tsv")
    results: dict[str, dict[str, str]] = {}
    for qid, (taxopath, conf) in assigned.items():
        # deepest taxopath element that the taxonomy knows -> the call
        nodes = [n for n in taxopath.split(";") if n]
        node = tlca.UNCLASSIFIED
        for n in reversed(nodes):
            if n in tlca.RETRO_PARENT:
                node = n
                break
        results[qid] = {
            "taxon_call": node,
            "rank": tlca.rank_of(node),
            "confidence": f"{conf:.3f}",
            "method": "placement",
        }
    return results
