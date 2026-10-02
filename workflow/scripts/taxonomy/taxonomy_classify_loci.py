"""Per-gene, mosaic-aware ERV locus classifier (trial).

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
import csv
import hashlib
import logging
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, NamedTuple

import pandas as pd
import taxonomy_lca as tlca
import taxonomy_placement
from Bio.Seq import Seq

from log import OK, job_logging, run_main
from tabular import gff3_attributes, gff3_features

# Domain-class semantics are shared with the scanner. The scanner imports THIS
# module for locus grouping, so the shared piece lives in its own module to keep
# that from becoming a cycle.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "domains"))
from domain_classes import load_classes
from domain_classes import summarise as summarise_domains

from external import run_tool

logger = logging.getLogger(__name__)

BLASTX = "blastx"  # translated search; all tools resolved from PATH (RetroSeek env)
MAKEBLASTDB = "makeblastdb"
# Per-locus marker regions are cut with a Bioconductor (Biostrings) helper rather
# than bedtools - RetroSeek keeps all range/sequence work inside Bioconductor.
_EXTRACT_R = Path(__file__).resolve().parent / "extract_region_fasta.R"

# The call a locus reports when none of its genes earned any call.
_UNCLASSIFIED_CALL = {
    "taxon_call": tlca.UNCLASSIFIED,
    "rank": "none",
    "confidence": "0.000",
    "method": "lca",
}
# The main gene set, the gene reliability order, the genomic gene order and the
# diagnostic genes are all derived at RUNTIME (from --main-probes, --gene-priority,
# --gene-order and the reference), never hard-coded, so the classifier is
# probe/gene-agnostic. See auto_diagnostic() and _assemble().


# ---------------------------------------------------------------- loci + regions
def parse_valid_full(gff3: Path) -> list[dict[str, str]]:
    """Read every per-hit feature of a valid-tier GFF3 as a flat string record.

    Rows that are not feature rows stop the job (`tabular.gff3_features`). The
    probe becomes the upper-cased ``gene`` (``OTHER`` when absent) and
    GFF3-escaped label separators are decoded. The strand test is a substring test against
    ``"+-"``: ``.`` and ``?`` read as ``+``, while an empty column (and the
    literal ``+-``) passes through unchanged.
    """
    return [_feature(fields) for fields, _, _ in gff3_features(gff3)]


def _feature(f: list[str]) -> dict[str, str]:
    """One valid-track feature from its nine GFF3 columns."""
    attributes = gff3_attributes(f[8])
    return {
        "seqname": f[0],
        "start": f[3],
        "end": f[4],
        "strand": f[6] if f[6] in "+-" else "+",
        "gene": attributes.get("probe", "OTHER").upper(),
        "parent": attributes.get("Parent", ""),
        "label": attributes.get("label", "").replace("%3b", ";"),
        # oversized rides the orphan track (overlap cluster wider than the
        # widest real provirus); the LTR-flanked track carries no attr.
        "oversized": attributes.get("oversized", "False"),
    }


def _group_by_parent(
    feats: list[dict[str, str]],
) -> dict[tuple[str, str], list[dict[str, str]]]:
    """Features keyed by (seqname, Parent); each parentless feature is its own group."""
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    orphan = 0
    for ft in feats:
        if ft["parent"]:
            key = (ft["seqname"], ft["parent"])
        else:
            key = (ft["seqname"], f"_orphan{orphan}")
            orphan += 1
        groups[key].append(ft)
    return groups


def _gene_spans(members: list[dict[str, str]]) -> dict[str, tuple[int, int]]:
    """Per gene, the span covering all of its hits in one locus (first-seen order)."""
    genes: dict[str, tuple[int, int]] = {}
    for m in members:
        s, e = int(m["start"]), int(m["end"])
        if m["gene"] in genes:
            old_s, old_e = genes[m["gene"]]
            s, e = min(old_s, s), max(old_e, e)
        genes[m["gene"]] = (s, e)
    return genes


def _strand_vote(members: list[dict[str, str]]) -> tuple[str, bool]:
    """The strand most hits of a locus are on, and whether that vote was tied."""
    votes = Counter(m["strand"] for m in members).most_common(2)
    tied = len(votes) == 2 and votes[0][1] == votes[1][1]
    return votes[0][0], tied


def _locus(
    index: int, seqname: str, parent: str, members: list[dict[str, str]]
) -> dict[str, Any]:
    """One locus from the features that share its Parent."""
    genes = _gene_spans(members)
    strand, strand_tie = _strand_vote(members)
    labels = {
        lab.strip() for m in members for lab in m["label"].split(";") if lab.strip()
    }
    return {
        "id": f"L{index}",
        "seqname": seqname,
        "parent": parent,
        "strand": strand,
        # A tied vote leaves the direction unknown: canonical_order stays blank.
        "strand_tie": strand_tie,
        "start": min(s for s, _ in genes.values()),
        "end": max(e for _, e in genes.values()),
        "genes": genes,
        # oversized: an orphan overlap-cluster wider than any real provirus;
        # members share a cluster, so any "True" marks the locus.
        "oversized": str(any(m.get("oversized", "False") == "True" for m in members)),
        "probe_label_set": ";".join(sorted(labels)),
    }


def build_loci(feats: list[dict[str, str]]) -> list[dict[str, Any]]:
    """Group features into LTR-element loci by their ``Parent=`` attribute; gene-partition each.

    Parent (the enclosing LTR_retrotransposon id) is emitted natively on the valid track by
    ``range_analysis/validation.R``. Parentless features become their own locus.
    """
    groups = _group_by_parent(feats)
    return [
        _locus(i, seqname, parent, members)
        for i, ((seqname, parent), members) in enumerate(groups.items())
    ]


def annotate_loci_with_domains(
    loci: list[dict[str, Any]],
    domains_parquet: Path,
    classes_tsv: Path,
    scanned_txt: Path,
) -> None:
    """Attach domain evidence to each locus in place, from the domain scan.

    Four columns are set: `domain_tier` (ADR-009 values, recomputed from curated
    classes rather than the retired name regexes), `domain_evidence` (the
    strongest class present), `domain_names` (the distinct families found) and
    `domain_source`.

    `domain_source` is the column that fixes the original defect. Before this,
    99.04% of `non_domain` rows meant "never assessed" yet were indistinguishable
    from the 292 that meant "assessed and empty". A locus present in the scanned
    roster but absent from the hit table is genuinely `non_domain`; one missing
    from the roster is `not_scanned`.

    Loci are joined on `{seqname}|{parent}`, the natural key. The catalog's own
    `L{i}` ids are positional and must never be used to join across processes.
    """
    hits = pd.read_parquet(domains_parquet).to_dict("records")
    per_locus = summarise_domains(hits, load_classes(classes_tsv))
    scanned = {
        line.strip()
        for line in scanned_txt.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    for lc in loci:
        key = f"{lc['seqname']}|{lc['parent']}"
        found = per_locus.get(key)
        if found is not None:
            lc.update(found)
            lc["domain_source"] = "scan"
        elif key in scanned:
            lc["domain_tier"] = "non_domain"
            lc["domain_evidence"] = "none"
            lc["domain_names"] = ""
            lc["domain_source"] = "scan"
        else:
            lc["domain_tier"] = "non_domain"
            lc["domain_evidence"] = "none"
            lc["domain_names"] = ""
            lc["domain_source"] = "not_scanned"


def write_region_bed(loci: list[dict[str, Any]], bed: Path) -> None:
    """Write one BED row per (locus, gene) region, named `{locus id}|{gene}`."""
    with bed.open("w", encoding="utf-8") as fh:
        for lc in loci:
            for gene, (s, e) in lc["genes"].items():
                fh.write(
                    f"{lc['seqname']}\t{max(0, s - 1)}\t{e}\t{lc['id']}|{gene}\t.\t{lc['strand']}\n"
                )


# ---------------------------------------------------------------- search
def build_db(faa: Path, db: Path) -> None:
    """Build the protein BLAST database from `faa`, unless it already exists."""
    if not db.with_suffix(".pin").exists():
        run_tool([MAKEBLASTDB, "-in", str(faa), "-dbtype", "prot", "-out", str(db)])


def search(query: Path, db: Path, out: Path, evalue: float, threads: int) -> None:
    """Run blastx of `query` against `db`, writing tabular hits to `out`."""
    run_tool(
        [
            BLASTX,
            "-query",
            str(query),
            "-db",
            str(db),
            "-out",
            str(out),
            "-outfmt",
            # qframe = query translation frame; pident = amino-acid identity of
            # the alignment, for the nearest-virus columns (ADR-024)
            "6 qseqid sseqid bitscore qframe pident",
            "-evalue",
            str(evalue),
            "-max_target_seqs",
            "25",
            "-num_threads",
            str(threads),
        ]
    )


def translate_frame(dna: str, frame: int) -> str:
    """Translate `dna` in blastx frame 1 to 3, or -1 to -3 (reverse strand)."""
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
    placement_out: Path | None = None,
    genome_name: str = "",
    domains_parquet: Path | None = None,
    domain_classes: Path | None = None,
    scanned_txt: Path | None = None,
    *,
    gene_priority: list[str] | None = None,
    gene_order: list[str] | None = None,
    placement_min_weight: float = 0.0,
) -> list[dict[str, str]]:
    """Classify every locus of one genome's valid track; one record per locus.

    Loci come from the track grouped by ``Parent=``. Each (locus, gene) region is
    cut from the genome and searched with blastx against the reference; the
    placement genes are also placed on their reference trees. The placement
    evidence is published to ``placement_out`` when given. Raises when an
    external tool fails (``run_tool``).
    """
    workdir.mkdir(parents=True, exist_ok=True)
    loci = build_loci(parse_valid_full(gff3))  # grouped by Parent= in the valid track
    if domains_parquet and domain_classes and scanned_txt:
        annotate_loci_with_domains(loci, domains_parquet, domain_classes, scanned_txt)
    fna, hits_path = _search_regions(loci, genome, ref_dir, workdir, evalue, threads)

    ref_csv = ref_dir / "retro_reference.csv"
    taxon_of, gene_of, axis = _ref_maps(ref_csv)
    diagnostic = auto_diagnostic(gene_of, taxon_of)  # genes in only one axis taxon
    region_seq = _load_regions(fna)
    virus_of = _virus_names(ref_csv)
    hits, best = _read_blastx(hits_path, taxon_of, virus_of)
    placement = _place_genes(
        placement_genes, best, region_seq, min_orf, ref_dir, workdir
    )
    if placement_out is not None:
        stem_prefix = f"{genome_name or genome.stem}.{source}"
        _export_placements(
            placement_genes, workdir, ref_dir, placement_out, stem_prefix
        )

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
        gene_priority=gene_priority,
        gene_order=gene_order,
        best=best,
        placement_min_weight=placement_min_weight,
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


def _virus_names(ref_csv: Path) -> dict[str, str]:
    """Each reference accession's virus: the last [bracketed] name of its defline.

    NCBI protein deflines end in the source organism, e.g. `gag protein [Mouse
    mammary tumor virus]`; an earlier bracket can hold a strain. A defline with
    no bracket is named by its accession, so the name is never blank.
    """
    names = {}
    with ref_csv.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            found = re.findall(r"\[([^\]]+)\]", r["defline"])
            names[r["accession"]] = found[-1] if found else r["accession"]
    return names


def auto_diagnostic(
    gene_of: dict[str, str], taxon_of: dict[str, str]
) -> dict[str, str]:
    """Genes whose reference members all belong to ONE axis taxon -> {gene: taxon} (presence-diagnostic).

    Data-derived (e.g. REX/TAX -> Deltaretrovirus), not hard-coded - so any probe set works.
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
    at the requested rank - the call itself when it is already there. Rank-
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


def _search_regions(
    loci: list[dict[str, Any]],
    genome: Path,
    ref_dir: Path,
    workdir: Path,
    evalue: float,
    threads: int,
) -> tuple[Path, Path]:
    """Cut each (locus, gene) region and blastx it; returns (regions FASTA, hits TSV)."""
    bed = workdir / "regions.bed"
    fna = workdir / "regions.fna"
    hits_path = workdir / "hits.tsv"
    db = workdir / "ref_db"
    write_region_bed(loci, bed)
    # strand-aware region extraction via Biostrings (Bioconductor), replacing
    # `bedtools getfasta -s -nameOnly`.
    run_tool(
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
    return fna, hits_path


class BestHit(NamedTuple):
    """A region's best blastx hit: its virus, identity, strength and frame."""

    virus: str
    identity: float  # percent amino-acid identity over the alignment
    bitscore: float
    frame: int  # blastx query frame, 1 to 3 or -1 to -3


def _read_blastx(
    hits_path: Path, taxon_of: dict[str, str], virus_of: dict[str, str]
) -> tuple[dict[str, list[tuple[str, float]]], dict[str, BestHit]]:
    """Per region (``locus|gene``): its (taxon, bitscore) hits, and its best hit.

    The best hit is the highest bitscore. It gives the frame placement translates
    in, and the nearest virus (ADR-024), named from ``virus_of``. Hits to
    accessions missing from the reference table are ignored.
    """
    hits: dict[str, list[tuple[str, float]]] = defaultdict(list)
    best: dict[str, BestHit] = {}
    with hits_path.open(encoding="utf-8") as fh:
        for line in fh:
            qid, sid, bits, frame, pident = line.rstrip("\n").split("\t")
            qid = qid.split("(")[0]
            accession = sid.split()[0]
            taxon = taxon_of.get(accession)
            if not taxon:
                continue
            b = float(bits)
            hits[qid].append((taxon, b))
            if qid not in best or b > best[qid].bitscore:
                # virus_of has every accession taxon_of has: one CSV feeds both.
                best[qid] = BestHit(virus_of[accession], float(pident), b, int(frame))
    return hits, best


def _placement_queries(
    gene: str,
    best: dict[str, BestHit],
    region_seq: dict[str, str],
    min_orf: int,
) -> dict[str, str]:
    """One gene's regions translated in their best blastx frame, stops as X.

    Translations shorter than ``min_orf`` residues are left out.
    """
    queries: dict[str, str] = {}
    for qid, hit in best.items():
        if qid.endswith(f"|{gene}") and qid in region_seq:
            prot = translate_frame(region_seq[qid], hit.frame).replace("*", "X")
            if len(prot) >= min_orf:
                queries[qid] = prot
    return queries


def _place_genes(
    placement_genes: set[str],
    best: dict[str, BestHit],
    region_seq: dict[str, str],
    min_orf: int,
    ref_dir: Path,
    workdir: Path,
) -> dict[str, dict[str, str]]:
    """Placement calls per region, one batch per placement gene."""
    placement: dict[str, dict[str, str]] = {}
    for gene in placement_genes:
        queries = _placement_queries(gene, best, region_seq, min_orf)
        if queries:
            placement.update(
                taxonomy_placement.place(
                    queries, ref_dir, gene, workdir / f"place_{gene}"
                )
            )
    return placement


def _export_placements(
    placement_genes: set[str],
    workdir: Path,
    ref_dir: Path,
    out_dir: Path,
    stem_prefix: str,
) -> None:
    """Publish the placement evidence before the scratch workdir is cleared.

    Runs for every placement gene, including those where `place()` never ran (no
    queries, all-gap alignment, missing tree package): export writes a valid
    empty jplace in that case so a rule declaring it still resolves.
    """
    for gene in sorted(placement_genes):
        taxonomy_placement.export_placement(
            workdir / f"place_{gene}", ref_dir, gene, out_dir, f"{stem_prefix}.{gene}"
        )


@dataclass(frozen=True)
class _Assembly:
    """The evidence and settings every locus record of one genome is built from."""

    hits: dict[str, list[tuple[str, float]]]
    placement: dict[str, dict[str, str]]
    diagnostic: dict[str, str]
    top_percent: float
    axis: set[str]
    main_probes: list[str]  # which genes are main; its length is the completeness count
    gene_priority: dict[str, int]  # position in the reliability list: lower wins
    gene_order: list[str]  # the genes 5' to 3', for canonical_order
    main_set: set[str]
    ref_version: str
    confidence_min: float
    structure_full_min: float
    source: str
    segment_rank: str
    best: dict[str, BestHit]  # region -> its best blastx hit (ADR-024)
    placement_min_weight: float  # below it a placement yields to blastx (ADR-025)


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
    *,
    gene_priority: list[str] | None = None,
    gene_order: list[str] | None = None,
    best: dict[str, BestHit] | None = None,
    placement_min_weight: float = 0.0,
) -> list[dict[str, str]]:
    """One output record per locus, from its per-gene evidence.

    The gene lists come from the user's config, never hard-coded (ADR-022):
    ``main_probes`` says which genes are main, ``gene_priority`` which gene's
    call wins (most reliable first), ``gene_order`` how the genes lie 5' to 3'.
    Either of the last two falls back to ``main_probes`` when not given.
    """
    asm = _Assembly(
        hits=hits,
        placement=placement,
        diagnostic=diagnostic,
        top_percent=top_percent,
        axis=axis,
        main_probes=main_probes,
        gene_priority={g: i for i, g in enumerate(gene_priority or main_probes)},
        gene_order=gene_order or main_probes,
        main_set=set(main_probes),
        ref_version=ref_version,
        confidence_min=confidence_min,
        structure_full_min=structure_full_min,
        source=source,
        segment_rank=segment_rank,
        best=best or {},
        placement_min_weight=placement_min_weight,
    )
    return [_locus_record(lc, asm) for lc in loci]


def _gene_call(qid: str, gene: str, asm: _Assembly) -> dict[str, str]:
    """The call for one gene region of a locus.

    Placement wins, but only when it resolves an axis taxon (ADR-008) with a
    placement weight (aLWR) of at least ``placement_min_weight`` (ADR-025); a
    diagnostic gene is called by presence; anything else by weighted-LCA over
    its blastx hits.
    """
    pl = asm.placement.get(qid)
    if (
        pl
        and pl["taxon_call"] in asm.axis
        and float(pl["confidence"]) >= asm.placement_min_weight
    ):
        return pl
    if gene in asm.diagnostic:
        node = asm.diagnostic[gene]
        return {
            "taxon_call": node,
            "rank": tlca.rank_of(node),
            "confidence": "1.000",
            "method": "presence",
        }
    node, wconf = tlca.weighted_lca(asm.hits.get(qid, []), asm.top_percent)
    return {
        "taxon_call": node,
        "rank": tlca.rank_of(node),
        "confidence": f"{wconf:.3f}",
        "method": "lca",
    }


def _priority(gene: str, asm: _Assembly) -> int:
    """The gene's place in ``gene_priority``; a gene it does not list comes last."""
    return asm.gene_priority.get(gene, len(asm.gene_priority))


def _call_gene(
    per_gene: dict[str, dict[str, str]],
    confident: dict[str, dict[str, str]],
    asm: _Assembly,
) -> str | None:
    """The gene whose call the locus reports, or None when no gene has one.

    Among genes resolved to an axis taxon: placement first, then marker
    reliability (``gene_priority``), then confidence. Raw confidence alone would
    favour ENV, because it depends on how much competes. With no gene resolved,
    the first gene with any call is used.
    """
    if confident:
        return min(
            confident,
            key=lambda g: (
                0 if confident[g]["method"] == "placement" else 1,
                _priority(g, asm),
                -float(confident[g]["confidence"]),
            ),
        )
    return next(
        (g for g, c in per_gene.items() if c["taxon_call"] != tlca.UNCLASSIFIED), None
    )


def _canonical_order(lc: dict[str, Any], gene_order: list[str]) -> str:
    """Whether the locus's genes lie in ``gene_order``, read along its strand.

    ``"True"`` or ``"False"``, or ``""`` when there is nothing to check: fewer
    than two of the listed genes are present, or the strand is unknown (its
    hits' vote was tied). Genes are ordered by their 5' end, which is the start
    coordinate on the plus strand and the end coordinate on the minus strand,
    so a provirus and its mirror image get the same answer even when one gene's
    span lies inside another's.
    """
    genes = lc["genes"]
    listed = [g for g in gene_order if g in genes]
    if len(listed) < 2 or lc.get("strand_tie") or lc["strand"] not in ("+", "-"):
        return ""
    if lc["strand"] == "+":
        along_strand = sorted(listed, key=lambda g: genes[g][0])
    else:
        along_strand = sorted(listed, key=lambda g: -genes[g][1])
    return str(along_strand == listed)


def _structure(lc: dict[str, Any], asm: _Assembly) -> dict[str, str]:
    """Structural columns over gene content: count, completeness, order, class.

    This loci table IS the genus-founded ERV assembly, so it carries the
    structure the legacy erv_like tier reported.
    """
    present_main = [g for g in asm.main_probes if g in lc["genes"]]
    completeness = len(present_main) / len(asm.main_probes) if asm.main_probes else 0.0
    # Discrete structural class over gene content (ADR-009): a single main gene
    # is a 'gene' fragment; a multi-gene locus is 'full' once its completeness
    # clears structure_full_min, else 'partial'. Deliberately gene-content only:
    # LTR-pair structure lives in the anchoring axis and the solo-LTR module.
    if len(present_main) <= 1:
        structure_class = "gene"
    elif completeness >= asm.structure_full_min:
        structure_class = "full"
    else:
        structure_class = "partial"
    return {
        "n_main_genes": str(len(present_main)),
        "completeness": f"{completeness:.3f}",
        "canonical_order": _canonical_order(lc, asm.gene_order),
        "structure_class": structure_class,
    }


NEAREST_COLUMNS = [
    "nearest_virus",
    "nearest_virus_identity",
    "nearest_virus_gene",
    "per_gene_nearest",
]


def _nearest(
    lc: dict[str, Any], asm: _Assembly, call_gene: str | None
) -> dict[str, str]:
    """The locus's nearest reference virus, with its amino-acid identity (ADR-024).

    Each gene region's best blastx hit names a reference virus. The headline gene
    is ``call_gene``, the gene behind the locus's taxon call, so the name and the
    call come from the same evidence. When that gene has no reference hit (or the
    locus has no call), ``gene_priority`` decides, then the strongest hit. The
    identity travels with the name, because a nearest virus is usually a distant
    relative, not an instance (measurements in ADR-024). All four columns are
    blank without a reference hit.
    """
    found = {g: hit for g in lc["genes"] if (hit := asm.best.get(f"{lc['id']}|{g}"))}
    if not found:
        return dict.fromkeys(NEAREST_COLUMNS, "")
    if call_gene is not None and call_gene in found:
        gene = call_gene
    else:
        gene = min(found, key=lambda g: (_priority(g, asm), -found[g].bitscore, g))
    return {
        "nearest_virus": found[gene].virus,
        "nearest_virus_identity": f"{found[gene].identity:.1f}",
        "nearest_virus_gene": gene,
        "per_gene_nearest": ";".join(
            f"{g}:{hit.virus}({hit.identity:.1f})" for g, hit in sorted(found.items())
        ),
    }


def _mosaic(
    confident: dict[str, dict[str, str]], main_set: set[str]
) -> tuple[str, str]:
    """(is_mosaic, composition): main genes resolved to more than one taxon.

    Only main genes count (OTHER is excluded; ENV is noisy but kept as a main
    gene); the composition lists every resolved gene.
    """
    distinct = {c["taxon_call"] for g, c in confident.items() if g in main_set}
    if len(distinct) <= 1:
        return "False", ""
    return "True", ";".join(
        f"{g}:{c['taxon_call']}" for g, c in sorted(confident.items())
    )


def _locus_record(lc: dict[str, Any], asm: _Assembly) -> dict[str, str]:
    """The output row for one locus; column order is the ``--out`` CSV's."""
    per_gene = {g: _gene_call(f"{lc['id']}|{g}", g, asm) for g in lc["genes"]}
    # 'confident' = resolved to an axis taxon (ADR-008), not the old rank=='genus'.
    confident = {g: c for g, c in per_gene.items() if c["taxon_call"] in asm.axis}
    call_gene = _call_gene(per_gene, confident, asm)
    call = per_gene[call_gene] if call_gene else dict(_UNCLASSIFIED_CALL)
    taxon_call = call["taxon_call"]
    is_mosaic, composition = _mosaic(confident, asm.main_set)
    # blastx evidence depth, summed over the locus's gene regions. Zero means valid
    # LTR structure but NO protein homology to the reference: the candidate-novel
    # retrovirus signal the loss analysis surfaces.
    n_blastx_hits = sum(len(asm.hits.get(f"{lc['id']}|{g}", [])) for g in lc["genes"])
    return {
        "id": lc["id"],
        "seqname": lc["seqname"],
        "start": str(lc["start"]),
        "end": str(lc["end"]),
        "strand": lc["strand"],
        "parent": lc["parent"],
        "genes_present": ",".join(sorted(lc["genes"])),
        **_structure(lc, asm),
        "domain_tier": lc.get("domain_tier", "non_domain"),
        "domain_evidence": lc.get("domain_evidence", "none"),
        "domain_names": lc.get("domain_names", ""),
        "domain_source": lc.get("domain_source", "not_scanned"),
        "oversized": lc.get("oversized", "False"),
        "taxon_call": taxon_call,
        "rank": call["rank"],
        # Rank roll-up for the segmentation stage (ADR-011): the locus's
        # ancestor at classification.segment_rank, or unassigned_at_<rank>.
        "segment": segment_of(taxon_call, asm.segment_rank),
        "segment_rank": asm.segment_rank,
        # resolved = the call landed on a declared axis taxon (ADR-008), vs an
        # honest LCA backoff to an interior ancestor.
        "resolved": str(taxon_call in asm.axis),
        "confidence": call["confidence"],
        # HC/LC against classification.confidence_min; the floor itself is HC.
        "confidence_tag": "LC"
        if float(call["confidence"]) < asm.confidence_min
        else "HC",
        "n_blastx_hits": str(n_blastx_hits),
        "method": call["method"],
        "per_gene": ";".join(
            f"{g}:{c['taxon_call']}({c['method']},{c['confidence']})"
            for g, c in sorted(per_gene.items())
        ),
        **_nearest(lc, asm, call_gene),
        "is_mosaic": is_mosaic,
        "mosaic_composition": composition,
        "erv_class": tlca.ERV_CLASS.get(taxon_call, ""),
        "probe_label_set": lc["probe_label_set"],
        "ref_version": asm.ref_version,
        "source": asm.source,
    }


def gate_classified(records: list[dict[str, str]]) -> list[dict[str, str]]:
    """Keep only records that earned a taxonomic call (drop UNCLASSIFIED).

    This is the orphan-recovery gate: a non-LTR-associated orphan is retained
    only if blastx resolved it to a taxon (axis member/backoff) - earning a
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


# Per-locus class counts: (metric, record column, value counted). Structural
# class and domain tier of the assembled loci (ADR-009); the grain is per locus,
# distinct from the per-hit valid_* tier counts emitted by ranges_analysis.R.
_CLASS_COUNTS = (
    ("structure_full", "structure_class", "full"),
    ("structure_partial", "structure_class", "partial"),
    ("structure_gene", "structure_class", "gene"),
    ("loci_domain_selected", "domain_tier", "domain_selected"),
    ("loci_domain_unlisted", "domain_tier", "domain_unlisted"),
    ("loci_non_domain", "domain_tier", "non_domain"),
)


def classification_counts(
    records: list[dict[str, str]], kept: list[dict[str, str]], source: str
) -> list[dict[str, Any]]:
    """Blastx-stage loss counters, as (metric, value) rows.

    They have the same shape as ranges_analysis.R's counts table, so the two
    UNION into one loss funnel.

    For the gated orphan run, ``records`` is the pre-gate set and ``kept`` the
    post-gate (recovered) set; for the ltr-flanked run the two are identical.
    """
    if source == "orphan":
        return [
            {"metric": "orphans_total", "value": len(records)},
            {"metric": "orphans_recovered", "value": len(kept)},
        ]
    classified = sum(r["taxon_call"] != tlca.UNCLASSIFIED for r in records)
    counts = [
        ("loci_total", len(records)),
        ("loci_classified", classified),
        ("loci_unclassified", len(records) - classified),
        ("loci_no_blastx_hit", sum(r["n_blastx_hits"] == "0" for r in records)),
    ]
    for metric, column, value in _CLASS_COUNTS:
        counts.append((metric, sum(r.get(column) == value for r in records)))
    return [{"metric": metric, "value": value} for metric, value in counts]


def summarise(records: list[dict[str, str]]) -> str:
    """A plain-text summary: counts, ranks, methods and resolved taxa of the calls."""
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
# Canonical schema for the loci table - fixed so an empty genome still writes a
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
    "domain_evidence",
    "domain_names",
    "domain_source",
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
    *NEAREST_COLUMNS,
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
                f"domain_evidence={r['domain_evidence']};"
                f"domain_source={r['domain_source']};"
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
        help="the main genes: the mosaic gene set and the completeness count. "
        "No gene is hard-coded.",
    )
    p.add_argument(
        "--gene-priority",
        default="",
        help="gene reliability for the locus call, most reliable first "
        "(default: the --main-probes order)",
    )
    p.add_argument(
        "--gene-order",
        default="",
        help="the genes 5' to 3', for canonical_order; any probe name "
        "(default: the --main-probes order)",
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
        "--placement-min-weight",
        type=float,
        default=0.0,
        help="placement weight (aLWR) a placement needs to win over the gene's "
        "blastx call; below it the blastx call stands. "
        "classification.placement_min_weight.",
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
        "--domains-parquet",
        type=Path,
        help="per-genome domain scan table (domains/scan_domains.py). Omit to "
        "leave every locus domain_source=not_scanned.",
    )
    p.add_argument(
        "--domain-classes",
        type=Path,
        help="curated Pfam class table, data/config/pfam_domain_classes.tsv.",
    )
    p.add_argument(
        "--domains-scanned",
        type=Path,
        help="roster of loci the scan covered; distinguishes non_domain from "
        "not_scanned.",
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
    # Ad-hoc single-CSV output (trial reproduction) - or the production output set:
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
    p.add_argument(
        "--out-placement-dir",
        type=Path,
        default=None,
        help=(
            "Directory to publish the placement evidence into: one "
            "{genome}.{tier}.{gene}.jplace + .labelled.newick per placement gene. "
            "Omit to leave the artifacts in the scratch workdir."
        ),
    )
    p.add_argument("--log", type=Path, help="job log (the Snakemake log: path)")
    return p


def _load_reference_taxonomy(ref_dir: Path) -> None:
    """Load the reference's own taxonomy and ERV classes when the build wrote them."""
    tax = ref_dir / "taxonomy.tsv"
    if tax.exists():
        tlca.load_taxonomy(tax)
    ervc = ref_dir / "erv_class.tsv"
    if ervc.exists():
        tlca.load_erv_class(ervc)


def _write_outputs(a: argparse.Namespace, records: list[dict[str, str]]) -> None:
    """Write whichever outputs were asked for: ad-hoc CSV, tables, IGV track."""
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
        logger.info("wrote track: %s", a.out_gff3)


def _comma_list(text: str) -> list[str]:
    """A comma-separated gene list from the command line, upper-cased, blanks dropped."""
    return [gene.strip().upper() for gene in text.split(",") if gene.strip()]


def _classify_from_args(a: argparse.Namespace) -> list[dict[str, str]]:
    """Run `classify` with the command line's inputs and settings."""
    return classify(
        a.gff3,
        a.genome,
        a.ref_dir,
        set(_comma_list(a.placement_genes)),
        _comma_list(a.main_probes),
        a.workdir,
        evalue=a.evalue,
        threads=a.threads,
        top_percent=a.top_percent,
        min_orf=a.min_orf,
        confidence_min=a.confidence_min,
        structure_full_min=a.structure_full_min,
        source=a.source,
        segment_rank=a.segment_rank,
        placement_out=a.out_placement_dir,
        genome_name=a.gff3.stem,
        domains_parquet=a.domains_parquet,
        domain_classes=a.domain_classes,
        scanned_txt=a.domains_scanned,
        gene_priority=_comma_list(a.gene_priority),
        gene_order=_comma_list(a.gene_order),
        placement_min_weight=a.placement_min_weight,
    )


def main() -> None:
    """Classify one genome's loci and write the requested outputs."""
    a = _build_arg_parser().parse_args()
    # One script, two rules (taxonomy_classify, taxonomy_orphans): the job log
    # path says which, and keeps parallel genomes apart.
    job_logging(a.log, f"taxonomy_classify_{a.source}")
    logger.info("classifying %s (source=%s)", a.gff3.stem, a.source)
    _load_reference_taxonomy(a.ref_dir)
    records = _classify_from_args(a)
    # Orphan-recovery gate: keep only loci that earned a taxonomic call. Counts
    # are computed over the PRE-gate set so the loss funnel can report what was
    # recovered vs. discarded. For the LTR-flanked run the gate is a no-op.
    kept = gate_classified(records) if a.gate_classified else records
    logger.info("classification summary\n%s", summarise(kept))
    if a.out_counts:
        write_counts(records, kept, a.source, a.out_counts)
        logger.info("wrote counts -> %s", a.out_counts)
    _write_outputs(a, kept)
    called = sum(1 for r in kept if r.get("taxon_call"))
    logger.log(OK, "%s loci, %s with a taxon call", f"{len(kept):,}", f"{called:,}")


if __name__ == "__main__":
    run_main(main)
