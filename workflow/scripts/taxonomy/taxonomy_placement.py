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

import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import taxonomy_lca as tlca
from Bio import SeqIO

MAFFT = "mafft"  # all tools resolved from PATH (the RetroSeek conda env)
EPA_NG = "epa-ng"
GAPPA = "gappa"

# The five per-placement values EPA-ng records, in jplace v3 order.
_JPLACE_FIELDS = [
    "edge_num",
    "likelihood",
    "like_weight_ratio",
    "distal_length",
    "pendant_length",
]


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


# ---------------------------------------------------------------------
# Artifact export
#
# EPA-ng and gappa write into a scratch workdir under data/tmp/, which is
# cleared between runs. The jplace is the evidence behind every taxon_call -
# which branch each locus attached to, with likelihood weights - and the
# standard interchange format iTOL and gappa read. These helpers copy it
# somewhere durable, and synthesise a valid empty file for the genomes where
# placement legitimately never ran, so a rule declaring it as an output still
# resolves. (tree_layout.py sets the same precedent with header-only CSVs.)
# ---------------------------------------------------------------------
def edge_numbered_newick(newick: str) -> str:
    """Return ``newick`` with a jplace ``{N}`` edge tag appended to every edge.

    jplace placements reference edges by number, so the embedded tree must carry
    them even when there are no placements to reference. Numbers are assigned in
    the order edges close, which is the post-order the format expects; with an
    empty placement list any self-consistent numbering parses correctly.
    """
    counter = 0

    def tag(_match: re.Match[str]) -> str:
        nonlocal counter
        out = f"{_match.group(0)}{{{counter}}}"
        counter += 1
        return out

    # An edge ends either at a branch length (":0.1") or at a bare node.
    body = newick.strip().rstrip(";")
    tagged = re.sub(r":-?[0-9.eE+-]+", tag, body)
    # The root edge carries no length of its own; give it the final number.
    return f"{tagged}{{{counter}}};"


def empty_jplace(newick: str, reason: str = "no queries were placed") -> str:
    """Return a valid jplace document with the real tree and no placements.

    The reference tree is carried through verbatim rather than replaced by a
    stand-in: downstream gappa commands read it to know what they are drawing,
    so a fabricated tree would silently mislabel a figure. ``reason`` is
    recorded because "zero placements" and "the stage never ran" are different
    claims and the file should say which one it is.
    """
    return json.dumps(
        {
            "version": 3,
            "fields": _JPLACE_FIELDS,
            "metadata": {"invocation": f"RetroSeek: {reason}"},
            "tree": edge_numbered_newick(newick),
            "placements": [],
        },
        indent=1,
    )


def _reference_newick(ref_dir: Path, gene: str) -> Path:
    """Path to the tree ``place()`` would actually run on, for this gene.

    Mirrors place()'s own preference for the raxml-ng-optimised tree, so the
    exported evidence matches the tree the calls were produced on.
    """
    trees = ref_dir / "trees"
    best = trees / f"{gene}.raxml.bestTree"
    fallback = trees / f"{gene}.treefile"
    if best.exists():
        return best
    if fallback.exists():
        return fallback
    raise FileNotFoundError(
        f"no reference tree for gene {gene!r} under {trees} "
        "(expected .raxml.bestTree or .treefile)"
    )


def export_placement(
    workdir: Path, ref_dir: Path, gene: str, out_dir: Path, stem: str
) -> tuple[Path, Path]:
    """Copy this gene's placement artifacts to ``out_dir`` as ``stem.*``.

    Returns ``(jplace_path, labelled_newick_path)``. When placement did not run -
    no queries carried the gene, every query aligned to all-gaps, or the gene had
    no tree package - an empty-but-valid jplace is written carrying the real
    reference tree, and the reference tree itself stands in for the labelled one.

    Raises ``FileNotFoundError`` when no reference tree exists at all, since
    there is then nothing honest to write.
    """
    reference = _reference_newick(ref_dir, gene)
    out_dir.mkdir(parents=True, exist_ok=True)
    jplace_out = out_dir / f"{stem}.jplace"
    newick_out = out_dir / f"{stem}.labelled.newick"

    jplace_src = workdir / "epa_result.jplace"
    newick_src = workdir / "labelled_tree.newick"

    if jplace_src.is_file():
        shutil.copyfile(jplace_src, jplace_out)
    else:
        jplace_out.write_text(
            empty_jplace(
                reference.read_text().strip(),
                reason=f"placement did not run for {gene}",
            )
        )
    shutil.copyfile(newick_src if newick_src.is_file() else reference, newick_out)
    return jplace_out, newick_out
