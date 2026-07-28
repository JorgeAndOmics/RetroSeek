"""
Per-gene, mosaic-aware ERV locus classifier (trial)
===================================================

Classifies ERV loci from their own sequence, per gene, with a method cascade:

    full valid GFF3 (per-hit features w/ probe=gene, Parent=LTR_retrotransposon)
      -> LTR-flanked loci (group by Parent; gene-partitioned regions)
      -> extract each (locus, gene) region (Biostrings; extract_region_fasta.R)
      -> blastx region vs independent reference  -> per-gene (taxon, bitscore) evidence
      -> per gene: placement (POL/GAG, if a tree exists) else weighted-LCA  [+ presence for REX/TAX]
      -> combine per-gene calls -> locus taxon_call + rank + confidence + method
         + is_mosaic + mosaic_composition + erv_class + detection provenance + ref_version

LCA is the universal default; placement is a dispatcher branch for configured genes.
Search defaults to blastx (no new dependency). See docs/taxonomy_classification/.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import logging
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
import taxonomy_lca as tlca
import taxonomy_placement
from Bio.Seq import Seq

from colored_logging import colored_logging

logger = logging.getLogger(__name__)

BLASTX = "blastx"  # translated search; all tools resolved from PATH (RetroSeek env)
MAKEBLASTDB = "makeblastdb"
# Per-locus marker regions are cut with a Bioconductor (Biostrings) helper rather
# than bedtools — RetroSeek keeps all range/sequence work inside Bioconductor.
_EXTRACT_R = Path(__file__).resolve().parent / "extract_region_fasta.R"

_PROBE = re.compile(r"probe=([^;\t]+)")
_PARENT = re.compile(r"Parent=([^;\t]+)")
_LABEL = re.compile(r"label=([^;\t]+)")
_ID = re.compile(r"ID=([^;\t]+)")
_DOMAIN_TIER = re.compile(r"domain_tier=([^;\t]+)")
_DOMAIN_HIT_CLASS = re.compile(r"domain_hit_class=([^;\t]+)")
_OVERSIZED = re.compile(r"oversized=([^;\t]+)")
# Per-provirus domain tier, strongest-wins when a locus's hits disagree (they
# shouldn't, since the tier is element-wise, but be defensive). See validation.R.
_DOMAIN_TIER_RANK = {"non_domain": 0, "domain_unlisted": 1, "domain_selected": 2}
# Gene reliability order, the mosaic gene set, and diagnostic genes are all derived at RUNTIME
# (from the user's ordered --main-probes and from the reference) — never hard-coded — so the
# classifier is probe/gene-agnostic. See auto_diagnostic() and _assemble().


def run(cmd: list[str], stdout: Any = None) -> None:
    res = subprocess.run(
        cmd,
        stdout=stdout or subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if res.returncode != 0:
        sys.stderr.write((res.stderr or "")[-2000:])
        raise SystemExit(f"command failed: {' '.join(cmd[:3])}…")


# ---------------------------------------------------------------- loci + regions
def parse_valid_full(gff3: Path) -> list[dict[str, str]]:
    feats: list[dict[str, str]] = []
    with gff3.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            probe = _PROBE.search(f[8])
            parent = _PARENT.search(f[8])
            label = _LABEL.search(f[8])
            tier = _DOMAIN_TIER.search(f[8])
            hit_class = _DOMAIN_HIT_CLASS.search(f[8])
            oversized = _OVERSIZED.search(f[8])
            feats.append(
                {
                    "seqname": f[0],
                    "start": f[3],
                    "end": f[4],
                    "strand": f[6] if f[6] in "+-" else "+",
                    "gene": (probe.group(1).upper() if probe else "OTHER"),
                    "parent": parent.group(1) if parent else "",
                    "label": (label.group(1) if label else "").replace("%3b", ";"),
                    # Domain labels ride the valid track (validation.R). The orphan
                    # track carries neither, so default to the weakest tier.
                    "domain_tier": tier.group(1) if tier else "non_domain",
                    "domain_hit_class": (
                        hit_class.group(1) if hit_class else "no_substring_match"
                    ),
                    # oversized rides the orphan track (overlap cluster wider than
                    # the widest real provirus); LTR-flanked track carries no attr.
                    "oversized": oversized.group(1) if oversized else "False",
                }
            )
    return feats


def load_elements(ltr_gff3: Path) -> dict[str, list[tuple[int, int, str]]]:
    """LTR_retrotransposon element ranges per seqname (sorted by start) from ltrdigest."""
    elems: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    with ltr_gff3.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "LTR_retrotransposon":
                continue
            eid = _ID.search(f[8])
            elems[f[0]].append(
                (int(f[3]), int(f[4]), eid.group(1) if eid else f"{f[0]}:{f[3]}")
            )
    for v in elems.values():
        v.sort()
    return elems


def assign_element(
    seqname: str,
    start: int,
    end: int,
    elems: dict[str, list[tuple[int, int, str]]],
    maxlen: dict[str, int],
) -> str:
    """Return the element id with greatest overlap of [start,end], or '' if none."""
    cand = elems.get(seqname)
    if not cand:
        return ""
    hi = bisect.bisect_right(
        [c[0] for c in cand], end
    )  # elements starting at/before end
    best, best_ov = "", 0
    lo = max(0, hi - 1)
    window_start = start - maxlen.get(seqname, 20000)
    while lo >= 0 and cand[lo][0] >= window_start:
        es, ee, eid = cand[lo]
        ov = min(end, ee) - max(start, es) + 1  # +1: GFF coords are 1-based inclusive
        if ov > best_ov:
            best, best_ov = eid, ov
        lo -= 1
    return best


def build_loci(feats: list[dict[str, str]]) -> list[dict[str, Any]]:
    """Group features into LTR-element loci by their ``Parent=`` attribute; gene-partition each.

    Parent (the enclosing LTR_retrotransposon id) is emitted natively on the valid track by
    ``range_analysis/validation.R``. Parentless features become their own locus.
    """
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    orphan = 0
    for ft in feats:
        if ft["parent"]:
            key = (ft["seqname"], ft["parent"])
        else:
            key = (ft["seqname"], f"_orphan{orphan}")
            orphan += 1
        groups[key].append(ft)

    loci: list[dict[str, Any]] = []
    for i, ((seqname, parent), members) in enumerate(groups.items()):
        strand = Counter(m["strand"] for m in members).most_common(1)[0][0]
        genes: dict[str, tuple[int, int]] = {}
        probe_labels: set[str] = set()
        # Per-provirus domain tier: element-wise, so a locus's members agree; take
        # the strongest defensively. All-orphan loci stay non_domain.
        domain_tier = "non_domain"
        # oversized: an orphan overlap-cluster wider than any real provirus; members
        # share a cluster, so any "True" marks the locus.
        oversized = "False"
        for m in members:
            if m.get("oversized", "False") == "True":
                oversized = "True"
            g = m["gene"]
            s, e = int(m["start"]), int(m["end"])
            if g in genes:
                genes[g] = (min(genes[g][0], s), max(genes[g][1], e))
            else:
                genes[g] = (s, e)
            for lab in m["label"].split(";"):
                if lab.strip():
                    probe_labels.add(lab.strip())
            mt = m.get("domain_tier", "non_domain")
            if _DOMAIN_TIER_RANK.get(mt, 0) > _DOMAIN_TIER_RANK.get(domain_tier, 0):
                domain_tier = mt
        start = min(v[0] for v in genes.values())
        end = max(v[1] for v in genes.values())
        loci.append(
            {
                "id": f"L{i}",
                "seqname": seqname,
                "parent": parent,
                "strand": strand,
                "start": start,
                "end": end,
                "genes": genes,
                "domain_tier": domain_tier,
                "oversized": oversized,
                "probe_label_set": ";".join(sorted(probe_labels)),
            }
        )
    return loci


def write_region_bed(loci: list[dict[str, Any]], bed: Path) -> None:
    with bed.open("w", encoding="utf-8") as fh:
        for lc in loci:
            for gene, (s, e) in lc["genes"].items():
                fh.write(
                    f"{lc['seqname']}\t{max(0, s - 1)}\t{e}\t{lc['id']}|{gene}\t.\t{lc['strand']}\n"
                )


# ---------------------------------------------------------------- search
def build_db(faa: Path, db: Path) -> None:
    if not db.with_suffix(".pin").exists():
        run([MAKEBLASTDB, "-in", str(faa), "-dbtype", "prot", "-out", str(db)])


def search(query: Path, db: Path, out: Path, evalue: float, threads: int) -> None:
    run(
        [
            BLASTX,
            "-query",
            str(query),
            "-db",
            str(db),
            "-out",
            str(out),
            "-outfmt",
            "6 qseqid sseqid bitscore qframe",  # qframe = query translation frame
            "-evalue",
            str(evalue),
            "-max_target_seqs",
            "25",
            "-num_threads",
            str(threads),
        ]
    )


def translate_frame(dna: str, frame: int) -> str:
    s = Seq(dna)
    if frame < 0:
        s = s.reverse_complement()  # type: ignore[no-untyped-call]
    s = s[abs(frame) - 1 :]
    s = s[: len(s) - (len(s) % 3)]
    return str(s.translate())  # type: ignore[no-untyped-call]


# ---------------------------------------------------------------- classify
def classify(
    gff3: Path,
    genome: Path,
    ref_dir: Path,
    placement_genes: set[str],
    main_probes: list[str],
    workdir: Path,
    evalue: float = 1e-3,
    threads: int = 1,
    top_percent: float = 0.10,
    min_orf: int = 30,
    confidence_min: float = 0.5,
    structure_full_min: float = 1.0,
    source: str = "ltr-flanked",
    segment_rank: str = "genus",
) -> list[dict[str, str]]:
    workdir.mkdir(parents=True, exist_ok=True)
    loci = build_loci(parse_valid_full(gff3))  # grouped by Parent= in the valid track
    bed = workdir / "regions.bed"
    fna = workdir / "regions.fna"
    db = workdir / "ref_db"
    hits_path = workdir / "hits.tsv"
    write_region_bed(loci, bed)
    # strand-aware region extraction via Biostrings (Bioconductor), replacing
    # `bedtools getfasta -s -nameOnly`.
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
    build_db(ref_dir / "retro_reference.faa", db)
    search(fna, db, hits_path, evalue, threads)

    taxon_of, gene_of, axis = _ref_maps(ref_dir / "retro_reference.csv")
    diagnostic = auto_diagnostic(
        gene_of, taxon_of
    )  # genes present in only one axis taxon
    region_seq = _load_regions(fna)
    # per region (locus|gene): taxon hits + best frame
    hits: dict[str, list[tuple[str, float]]] = defaultdict(list)
    best_frame: dict[str, tuple[float, int]] = {}
    with hits_path.open(encoding="utf-8") as fh:
        for line in fh:
            qid, sid, bits, frame = line.rstrip("\n").split("\t")
            qid = qid.split("(")[0]
            taxon = taxon_of.get(sid.split()[0])
            if not taxon:
                continue
            b = float(bits)
            hits[qid].append((taxon, b))
            if qid not in best_frame or b > best_frame[qid][0]:
                best_frame[qid] = (b, int(frame))

    # placement: batch per placement-gene
    placement: dict[str, dict[str, str]] = {}
    for gene in placement_genes:
        queries: dict[str, str] = {}
        for qid, (_, fr) in best_frame.items():
            if qid.endswith(f"|{gene}") and qid in region_seq:
                prot = translate_frame(region_seq[qid], fr).replace("*", "X")
                if len(prot) >= min_orf:
                    queries[qid] = prot
        if queries:
            res = taxonomy_placement.place(
                queries, ref_dir, gene, workdir / f"place_{gene}"
            )
            placement.update(res)

    return _assemble(
        loci,
        hits,
        placement,
        _ref_version(ref_dir),
        main_probes,
        diagnostic,
        top_percent,
        axis,
        confidence_min=confidence_min,
        structure_full_min=structure_full_min,
        source=source,
        segment_rank=segment_rank,
    )


def _ref_maps(ref_csv: Path) -> tuple[dict[str, str], dict[str, str], set[str]]:
    """Return (accession->taxon, accession->gene, axis) from the reference CSV.

    ``axis`` is the set of declared reference taxa (any rank); a locus is 'resolved'
    when its call lands on an axis member (ADR-008), replacing the old rank=='genus' test.
    """
    taxon_of, gene_of = {}, {}
    with ref_csv.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            taxon_of[r["accession"]] = r["taxon"]
            gene_of[r["accession"]] = r["gene"]
    return taxon_of, gene_of, set(taxon_of.values())


def auto_diagnostic(
    gene_of: dict[str, str], taxon_of: dict[str, str]
) -> dict[str, str]:
    """Genes whose reference members all belong to ONE axis taxon -> {gene: taxon} (presence-diagnostic).

    Data-derived (e.g. REX/TAX -> Deltaretrovirus), not hard-coded — so any probe set works.
    """
    by_gene: dict[str, set[str]] = defaultdict(set)
    for acc, gene in gene_of.items():
        by_gene[gene].add(taxon_of[acc])
    return {
        g: next(iter(gs)) for g, gs in by_gene.items() if len(gs) == 1 and g != "OTHER"
    }


def _load_regions(fna: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur: str | None = None
    buf: list[str] = []
    with fna.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buf)
                cur = line[1:].strip().split("(")[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur:
        seqs[cur] = "".join(buf)
    return seqs


def _ref_version(ref_dir: Path) -> str:
    h = hashlib.md5()
    h.update((ref_dir / "retro_reference.csv").read_bytes())
    return h.hexdigest()[:12]


def segment_of(taxon_call: str, segment_rank: str) -> str:
    """Roll a taxon call up to ``segment_rank`` (ADR-011).

    Walks the taxonomy from the call toward the root and returns the first node
    at the requested rank — the call itself when it is already there. Rank-
    agnostic by construction: ``segment_rank`` is any NCBI rank string, and no
    taxon name is ever hard-coded, so segmenting by family works exactly like
    segmenting by genus.

    A call ABOVE the requested rank (e.g. ``Retroviridae`` when segmenting by
    genus) has no genus ancestor and cannot be given one: it returns
    ``unassigned_at_<rank>`` rather than inventing precision the evidence does
    not support.
    """
    if not segment_rank:
        return ""
    for node in tlca.ancestors(taxon_call):
        if tlca.rank_of(node) == segment_rank:
            return node
    return f"unassigned_at_{segment_rank}"


def _assemble(
    loci: list[dict[str, Any]],
    hits: dict[str, list[tuple[str, float]]],
    placement: dict[str, dict[str, str]],
    ref_version: str,
    main_probes: list[str],
    diagnostic: dict[str, str],
    top_percent: float,
    axis: set[str],
    confidence_min: float = 0.5,
    structure_full_min: float = 1.0,
    source: str = "ltr-flanked",
    segment_rank: str = "genus",
) -> list[dict[str, str]]:
    # gene reliability + mosaic set derived from the user's ordered main_probes (no hard-coding)
    gene_priority = {g: i for i, g in enumerate(main_probes)}
    default_priority = len(main_probes) + 1
    main_set = set(main_probes)
    records: list[dict[str, str]] = []
    for lc in loci:
        per_gene: dict[str, dict[str, str]] = {}
        for gene in lc["genes"]:
            qid = f"{lc['id']}|{gene}"
            pl = placement.get(qid)
            if (
                pl and pl["taxon_call"] in axis
            ):  # placement refines ONLY when it resolves an axis taxon (ADR-008)
                per_gene[gene] = pl
            elif gene in diagnostic:
                node = diagnostic[gene]
                per_gene[gene] = {
                    "taxon_call": node,
                    "rank": tlca.rank_of(node),
                    "confidence": "1.000",
                    "method": "presence",
                }
            else:
                node, wconf = tlca.weighted_lca(hits.get(qid, []), top_percent)
                per_gene[gene] = {
                    "taxon_call": node,
                    "rank": tlca.rank_of(node),
                    "confidence": f"{wconf:.3f}",
                    "method": "lca",
                }
        # locus summary: choose by marker reliability (POL>GAG>…>ENV), placement preferred,
        # then confidence — NOT raw confidence (which is competition-dependent and favours ENV).
        # 'confident' = resolved to an axis taxon (ADR-008), replacing the old rank=='genus'.
        confident = {g: c for g, c in per_gene.items() if c["taxon_call"] in axis}
        if confident:
            best_gene = min(
                confident,
                key=lambda g: (
                    0 if confident[g]["method"] == "placement" else 1,
                    gene_priority.get(g, default_priority),
                    -float(confident[g]["confidence"]),
                ),
            )
            call = confident[best_gene]
            taxon_call, rank, conf, method = (
                call["taxon_call"],
                call["rank"],
                call["confidence"],
                call["method"],
            )
        else:
            # no axis taxon resolved -> best non-unclassified by rank, else unclassified
            ranked = [
                c for c in per_gene.values() if c["taxon_call"] != tlca.UNCLASSIFIED
            ]
            call = (
                ranked[0]
                if ranked
                else {
                    "taxon_call": tlca.UNCLASSIFIED,
                    "rank": "none",
                    "confidence": "0.000",
                    "method": "lca",
                }
            )
            taxon_call, rank, conf, method = (
                call["taxon_call"],
                call["rank"],
                call["confidence"],
                call["method"],
            )
        # mosaic only over main genes (exclude OTHER; ENV noisy but kept as a main gene)
        distinct = {c["taxon_call"] for g, c in confident.items() if g in main_set}
        # structural metrics — this loci table IS the genus-founded ERV assembly, so it
        # carries the same structure the legacy erv_like tier reported: how many main
        # genes are present (completeness) and whether they sit in canonical genomic order.
        present_main = [g for g in main_probes if g in lc["genes"]]
        by_pos = [
            g
            for g, _ in sorted(
                ((g, lc["genes"][g][0]) for g in present_main), key=lambda x: x[1]
            )
        ]
        canonical = bool(present_main) and by_pos in (present_main, present_main[::-1])
        # Discrete structural class over gene content (ADR-009): a single main
        # gene is a 'gene' fragment; a multi-gene locus is 'full' once its
        # completeness clears structure_full_min, else 'partial'. Deliberately
        # gene-content only — LTR-pair structure lives in the anchoring axis and
        # the solo-LTR module, not here.
        completeness_val = len(present_main) / len(main_probes) if main_probes else 0.0
        if len(present_main) <= 1:
            structure_class = "gene"
        elif completeness_val >= structure_full_min:
            structure_class = "full"
        else:
            structure_class = "partial"
        # blastx evidence depth for this locus (summed over its gene regions). Zero
        # means the locus carries valid LTR structure but NO protein homology to the
        # reference — the candidate-novel-retrovirus signal the loss analysis surfaces.
        n_blastx_hits = sum(len(hits.get(f"{lc['id']}|{g}", [])) for g in lc["genes"])
        # confidence tag: HC/LC against a user-adjustable floor (classification.confidence_min).
        # Threshold is inclusive — conf == floor is still High Confidence.
        confidence_tag = "LC" if float(conf) < confidence_min else "HC"
        records.append(
            {
                "id": lc["id"],
                "seqname": lc["seqname"],
                "start": str(lc["start"]),
                "end": str(lc["end"]),
                "strand": lc["strand"],
                "parent": lc["parent"],
                "genes_present": ",".join(sorted(lc["genes"])),
                "n_main_genes": str(len(present_main)),
                "completeness": f"{completeness_val:.3f}",
                "canonical_order": str(canonical),
                "structure_class": structure_class,
                "domain_tier": lc.get("domain_tier", "non_domain"),
                "oversized": lc.get("oversized", "False"),
                "taxon_call": taxon_call,
                "rank": rank,
                # Rank roll-up for the segmentation stage (ADR-011): the locus's
                # ancestor at classification.segment_rank, or unassigned_at_<rank>
                # when the call is coarser than that rank.
                "segment": segment_of(taxon_call, segment_rank),
                "segment_rank": segment_rank,
                # resolved = the call landed on a declared axis taxon (ADR-008),
                # vs an honest LCA-backoff to an interior ancestor. Rank-agnostic
                # 'confident' flag for downstream consumers (default axis=genera
                # -> resolved iff rank=='genus', so plots stay identical).
                "resolved": str(taxon_call in axis),
                "confidence": conf,
                "confidence_tag": confidence_tag,
                "n_blastx_hits": str(n_blastx_hits),
                "method": method,
                "per_gene": ";".join(
                    f"{g}:{c['taxon_call']}({c['method']},{c['confidence']})"
                    for g, c in sorted(per_gene.items())
                ),
                "is_mosaic": str(len(distinct) > 1),
                "mosaic_composition": (
                    ";".join(
                        f"{g}:{c['taxon_call']}" for g, c in sorted(confident.items())
                    )
                    if len(distinct) > 1
                    else ""
                ),
                "erv_class": tlca.ERV_CLASS.get(taxon_call, ""),
                "probe_label_set": lc["probe_label_set"],
                "ref_version": ref_version,
                "source": source,
            }
        )
    return records


def gate_classified(records: list[dict[str, str]]) -> list[dict[str, str]]:
    """Keep only records that earned a taxonomic call (drop UNCLASSIFIED).

    This is the orphan-recovery gate: a non-LTR-associated orphan is retained
    only if blastx resolved it to a taxon (axis member/backoff) — earning a
    classification *is* the evidence it is a real (possibly novel) retroviral orphan.
    """
    return [r for r in records if r["taxon_call"] != tlca.UNCLASSIFIED]


def write_counts(
    records: list[dict[str, str]],
    kept: list[dict[str, str]],
    source: str,
    out_counts: Path,
) -> None:
    """Write the blastx-stage loss counts CSV (metric,value)."""
    out_counts.parent.mkdir(parents=True, exist_ok=True)
    with out_counts.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["metric", "value"])
        w.writeheader()
        w.writerows(classification_counts(records, kept, source))


def classification_counts(
    records: list[dict[str, str]], kept: list[dict[str, str]], source: str
) -> list[dict[str, Any]]:
    """Blastx-stage loss counters, emitted in the same (metric, value) shape as
    ranges_analysis.R's counts table so the two UNION into one loss funnel.

    For the gated orphan run, ``records`` is the pre-gate set and ``kept`` the
    post-gate (recovered) set; for the ltr-flanked run the two are identical.
    """
    if source == "orphan":
        return [
            {"metric": "orphans_total", "value": len(records)},
            {"metric": "orphans_recovered", "value": len(kept)},
        ]
    return [
        {"metric": "loci_total", "value": len(records)},
        {
            "metric": "loci_classified",
            "value": sum(r["taxon_call"] != tlca.UNCLASSIFIED for r in records),
        },
        {
            "metric": "loci_unclassified",
            "value": sum(r["taxon_call"] == tlca.UNCLASSIFIED for r in records),
        },
        {
            "metric": "loci_no_blastx_hit",
            "value": sum(r["n_blastx_hits"] == "0" for r in records),
        },
        # Per-provirus structural class + domain tier of the assembled loci
        # (ADR-009). Grain is per-locus, distinct from the per-hit valid_* tier
        # counts emitted by ranges_analysis.R.
        {
            "metric": "structure_full",
            "value": sum(r.get("structure_class") == "full" for r in records),
        },
        {
            "metric": "structure_partial",
            "value": sum(r.get("structure_class") == "partial" for r in records),
        },
        {
            "metric": "structure_gene",
            "value": sum(r.get("structure_class") == "gene" for r in records),
        },
        {
            "metric": "loci_domain_selected",
            "value": sum(r.get("domain_tier") == "domain_selected" for r in records),
        },
        {
            "metric": "loci_domain_unlisted",
            "value": sum(r.get("domain_tier") == "domain_unlisted" for r in records),
        },
        {
            "metric": "loci_non_domain",
            "value": sum(r.get("domain_tier") == "non_domain" for r in records),
        },
    ]


def summarise(records: list[dict[str, str]]) -> str:
    total = len(records)
    placed = [r for r in records if r["taxon_call"] != tlca.UNCLASSIFIED]
    by_rank = Counter(r["rank"] for r in placed)
    taxa = Counter(r["taxon_call"] for r in placed)
    methods = Counter(r["method"] for r in placed)
    mosaic = sum(r["is_mosaic"] == "True" for r in records)
    lines = [
        f"loci: {total}   classified: {len(placed)}   mosaic: {mosaic}",
        "rank: " + ", ".join(f"{k}={v}" for k, v in by_rank.items()),
        "call method: " + ", ".join(f"{k}={v}" for k, v in methods.items()),
        "resolved taxa:",
    ]
    for t, n in taxa.most_common():
        lines.append(f"  {t:18s} {n:6d}  [{tlca.ERV_CLASS.get(t, '')}]")
    return "\n".join(lines)


# ---------------------------------------------------------------- outputs
# Canonical schema for the loci table — fixed so an empty genome still writes a
# well-formed parquet/csv (Snakemake output contract) instead of a headerless file.
LOCI_COLUMNS = [
    "id",
    "seqname",
    "start",
    "end",
    "strand",
    "parent",
    "genes_present",
    "n_main_genes",
    "completeness",
    "canonical_order",
    "structure_class",
    "domain_tier",
    "oversized",
    "taxon_call",
    "rank",
    "segment",
    "segment_rank",
    "resolved",
    "confidence",
    "confidence_tag",
    "n_blastx_hits",
    "method",
    "per_gene",
    "is_mosaic",
    "mosaic_composition",
    "erv_class",
    "probe_label_set",
    "ref_version",
    "source",
]


def write_tables(records: list[dict[str, str]], parquet: Path, csv_path: Path) -> None:
    """Write the per-locus taxon-call table as dual parquet + csv (always, even if empty)."""
    df = pd.DataFrame(records, columns=LOCI_COLUMNS)
    parquet.parent.mkdir(parents=True, exist_ok=True)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(parquet, index=False)
    df.to_csv(csv_path, index=False)


def write_track(records: list[dict[str, str]], gff3: Path, bed: Path) -> None:
    """Project per-locus calls to genome coordinates: IGV GFF3 + BED (colour-by-taxon)."""
    gff3.parent.mkdir(parents=True, exist_ok=True)
    with gff3.open("w", encoding="utf-8") as g, bed.open("w", encoding="utf-8") as b:
        g.write("##gff-version 3\n")
        for r in records:
            attrs = (
                f"ID={r['id']};taxon={r['taxon_call']};rank={r['rank']};"
                f"method={r['method']};confidence={r['confidence']};"
                f"confidence_tag={r['confidence_tag']};"
                f"structure_class={r['structure_class']};"
                f"domain_tier={r['domain_tier']};"
                f"mosaic={r['is_mosaic']};erv_class={r['erv_class']};"
                f"genes={r['genes_present']}"
            )
            g.write(
                f"{r['seqname']}\tRetroSeek\tERV_locus\t{r['start']}\t{r['end']}\t"
                f"{r['confidence']}\t{r['strand']}\t.\t{attrs}\n"
            )
            # BED is 0-based half-open; name carries the taxon call for IGV colour-by-name.
            b.write(
                f"{r['seqname']}\t{int(r['start']) - 1}\t{r['end']}\t"
                f"{r['id']}|{r['taxon_call']}\t0\t{r['strand']}\n"
            )


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Per-gene mosaic-aware ERV classifier")
    p.add_argument(
        "gff3", type=Path, help="FULL valid-tier GFF3 (per-hit, with Parent)"
    )
    p.add_argument("genome", type=Path, help="genome FASTA")
    p.add_argument("--ref-dir", type=Path, default=Path("data/taxonomy_reference"))
    p.add_argument(
        "--main-probes",
        default="POL,GAG,ENV,PRO",
        help="ordered gene reliability for the locus call + the mosaic gene set "
        "(earlier = preferred). Drives gene-agnostic behaviour; no gene hard-coded.",
    )
    p.add_argument(
        "--placement-genes",
        default="POL",
        help="genes classified by phylogenetic placement (others: weighted-LCA)",
    )
    p.add_argument("--evalue", type=float, default=1e-3, help="blastx e-value cutoff")
    p.add_argument(
        "--top-percent",
        type=float,
        default=0.10,
        help="weighted-LCA: keep hits within this fraction of the best bitscore",
    )
    p.add_argument(
        "--min-orf",
        type=int,
        default=30,
        help="min translated marker length (aa) to place",
    )
    p.add_argument(
        "--confidence-min",
        type=float,
        default=0.5,
        help="confidence floor below which a locus call is tagged 'LC' (low "
        "confidence); at or above it is 'HC'. classification.confidence_min.",
    )
    p.add_argument(
        "--structure-full-min",
        type=float,
        default=1.0,
        help="min gene completeness (fraction of main genes) for structure_class "
        "'full'; a single-gene locus is 'gene', below-threshold multi-gene is "
        "'partial'. classification.structure_full_min.",
    )
    p.add_argument(
        "--segment-rank",
        default="genus",
        help="taxonomic rank the `segment` column rolls each call up to "
        "(any NCBI rank: genus, subfamily, family...). Calls coarser than this "
        "rank become unassigned_at_<rank>. classification.segment_rank.",
    )
    p.add_argument(
        "--source",
        default="ltr-flanked",
        help="provenance stamp for every record ('ltr-flanked' LTR loci vs "
        "recovered 'orphan'). Lets downstream union/report split the two tiers.",
    )
    p.add_argument(
        "--gate-classified",
        action="store_true",
        help="drop UNCLASSIFIED records before writing (the orphan-recovery "
        "gate: keep an orphan only if it earned a taxonomic call).",
    )
    p.add_argument("--threads", type=int, default=1, help="blastx threads")
    p.add_argument("--workdir", type=Path, default=Path("/tmp/taxonomy_classify"))
    # Ad-hoc single-CSV output (trial reproduction) — or the production output set:
    p.add_argument("--out", type=Path, default=None, help="ad-hoc single CSV output")
    p.add_argument("--out-parquet", type=Path, default=None, help="loci table parquet")
    p.add_argument("--out-csv", type=Path, default=None, help="loci table csv")
    p.add_argument(
        "--out-gff3", type=Path, default=None, help="genome-coordinate IGV track"
    )
    p.add_argument(
        "--out-bed", type=Path, default=None, help="genome-coordinate BED track"
    )
    p.add_argument(
        "--out-counts",
        type=Path,
        default=None,
        help="blastx-stage loss counts CSV (metric,value) for the loss funnel",
    )
    return p


def main() -> int:
    a = _build_arg_parser().parse_args()

    # Per-(genome, tier) log file so the parallel per-genome invocations don't
    # clobber one another's log.
    colored_logging(log_file_name=f"taxonomy_classify_{a.gff3.stem}_{a.source}.txt")
    logger.info("classifying %s (source=%s)", a.gff3.stem, a.source)

    tax = a.ref_dir / "taxonomy.tsv"
    if tax.exists():
        tlca.load_taxonomy(tax)
    ervc = a.ref_dir / "erv_class.tsv"
    if ervc.exists():
        tlca.load_erv_class(ervc)
    pgenes = {g.strip().upper() for g in a.placement_genes.split(",") if g.strip()}
    main_probes = [g.strip().upper() for g in a.main_probes.split(",") if g.strip()]

    records = classify(
        a.gff3,
        a.genome,
        a.ref_dir,
        pgenes,
        main_probes,
        a.workdir,
        evalue=a.evalue,
        threads=a.threads,
        top_percent=a.top_percent,
        min_orf=a.min_orf,
        confidence_min=a.confidence_min,
        structure_full_min=a.structure_full_min,
        source=a.source,
        segment_rank=a.segment_rank,
    )
    # Orphan-recovery gate: keep only loci that earned a taxonomic call. Counts
    # are computed over the PRE-gate set so the loss funnel can report what was
    # recovered vs. discarded. For the LTR-flanked run the gate is a no-op.
    kept = gate_classified(records) if a.gate_classified else records
    logger.info("classification summary\n%s", summarise(kept))
    if a.out_counts:
        write_counts(records, kept, a.source, a.out_counts)
        logger.info("wrote counts -> %s", a.out_counts)
    records = kept

    if a.out and records:  # ad-hoc single CSV (back-compat for trial docs)
        with a.out.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(records[0].keys()))
            w.writeheader()
            w.writerows(records)
        logger.info("wrote -> %s", a.out)
    if a.out_parquet and a.out_csv:  # production dual-table output
        write_tables(records, a.out_parquet, a.out_csv)
        logger.info("wrote tables -> %s, %s", a.out_parquet, a.out_csv)
    if a.out_gff3 and a.out_bed:  # production IGV track
        write_track(records, a.out_gff3, a.out_bed)
        logger.info("wrote track -> %s", a.out_gff3)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
