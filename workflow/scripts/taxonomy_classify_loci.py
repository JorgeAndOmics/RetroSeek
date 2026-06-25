"""
Per-gene, mosaic-aware ERV locus classifier (trial)
===================================================

Classifies ERV loci from their own sequence, per gene, with a method cascade:

    full valid GFF3 (per-hit features w/ probe=gene, Parent=LTR_retrotransposon)
      -> LTR-element-anchored loci (group by Parent; gene-partitioned regions)
      -> extract each (locus, gene) region (bedtools getfasta)
      -> blastx region vs independent reference  -> per-gene (genus, bitscore) evidence
      -> per gene: placement (POL/GAG, if a tree exists) else weighted-LCA  [+ presence for REX/TAX]
      -> combine per-gene calls -> locus genus_call + rank + confidence + method
         + is_mosaic + mosaic_composition + erv_class + detection provenance + ref_version

LCA is the universal default; placement is a dispatcher branch for configured genes.
Search defaults to blastx (no new dependency). See docs/taxonomy_classification/.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
from Bio.Seq import Seq

import taxonomy_lca as tlca
import taxonomy_placement

BLASTX = "blastx"  # translated search; all tools resolved from PATH (RetroSeek env)
MAKEBLASTDB = "makeblastdb"
BEDTOOLS = "bedtools"

_PROBE = re.compile(r"probe=([^;\t]+)")
_PARENT = re.compile(r"Parent=([^;\t]+)")
_LABEL = re.compile(r"label=([^;\t]+)")
_ID = re.compile(r"ID=([^;\t]+)")
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
            feats.append(
                {
                    "seqname": f[0],
                    "start": f[3],
                    "end": f[4],
                    "strand": f[6] if f[6] in "+-" else "+",
                    "gene": (probe.group(1).upper() if probe else "OTHER"),
                    "parent": parent.group(1) if parent else "",
                    "label": (label.group(1) if label else "").replace("%3b", ";"),
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
        for m in members:
            g = m["gene"]
            s, e = int(m["start"]), int(m["end"])
            if g in genes:
                genes[g] = (min(genes[g][0], s), max(genes[g][1], e))
            else:
                genes[g] = (s, e)
            for lab in m["label"].split(";"):
                if lab.strip():
                    probe_labels.add(lab.strip())
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
) -> list[dict[str, str]]:
    workdir.mkdir(parents=True, exist_ok=True)
    loci = build_loci(parse_valid_full(gff3))  # grouped by Parent= in the valid track
    bed = workdir / "regions.bed"
    fna = workdir / "regions.fna"
    db = workdir / "ref_db"
    hits_path = workdir / "hits.tsv"
    write_region_bed(loci, bed)
    run(
        [
            BEDTOOLS,
            "getfasta",
            "-fi",
            str(genome),
            "-bed",
            str(bed),
            "-s",
            "-nameOnly",
            "-fo",
            str(fna),
        ]
    )
    build_db(ref_dir / "retro_reference.faa", db)
    search(fna, db, hits_path, evalue, threads)

    genus_of, gene_of = _ref_maps(ref_dir / "retro_reference.csv")
    diagnostic = auto_diagnostic(gene_of, genus_of)  # genes present in only one genus
    region_seq = _load_regions(fna)
    # per region (locus|gene): genus hits + best frame
    hits: dict[str, list[tuple[str, float]]] = defaultdict(list)
    best_frame: dict[str, tuple[float, int]] = {}
    with hits_path.open(encoding="utf-8") as fh:
        for line in fh:
            qid, sid, bits, frame = line.rstrip("\n").split("\t")
            qid = qid.split("(")[0]
            genus = genus_of.get(sid.split()[0])
            if not genus:
                continue
            b = float(bits)
            hits[qid].append((genus, b))
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
    )


def _ref_maps(ref_csv: Path) -> tuple[dict[str, str], dict[str, str]]:
    genus_of, gene_of = {}, {}
    with ref_csv.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            genus_of[r["accession"]] = r["genus"]
            gene_of[r["accession"]] = r["gene"]
    return genus_of, gene_of


def auto_diagnostic(
    gene_of: dict[str, str], genus_of: dict[str, str]
) -> dict[str, str]:
    """Genes whose reference members all belong to ONE genus -> {gene: genus} (presence-diagnostic).

    Data-derived (e.g. REX/TAX -> Deltaretrovirus), not hard-coded — so any probe set works.
    """
    by_gene: dict[str, set[str]] = defaultdict(set)
    for acc, gene in gene_of.items():
        by_gene[gene].add(genus_of[acc])
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


def _assemble(
    loci: list[dict[str, Any]],
    hits: dict[str, list[tuple[str, float]]],
    placement: dict[str, dict[str, str]],
    ref_version: str,
    main_probes: list[str],
    diagnostic: dict[str, str],
    top_percent: float,
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
                pl and pl["rank"] == "genus"
            ):  # placement refines ONLY when it resolves a genus
                per_gene[gene] = pl
            elif gene in diagnostic:
                node = diagnostic[gene]
                per_gene[gene] = {
                    "genus_call": node,
                    "rank": tlca.rank_of(node),
                    "confidence": "1.000",
                    "method": "presence",
                }
            else:
                node, wconf = tlca.weighted_lca(hits.get(qid, []), top_percent)
                per_gene[gene] = {
                    "genus_call": node,
                    "rank": tlca.rank_of(node),
                    "confidence": f"{wconf:.3f}",
                    "method": "lca",
                }
        # locus summary: choose by marker reliability (POL>GAG>…>ENV), placement preferred,
        # then confidence — NOT raw confidence (which is competition-dependent and favours ENV).
        confident = {g: c for g, c in per_gene.items() if c["rank"] == "genus"}
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
            genus_call, rank, conf, method = (
                call["genus_call"],
                call["rank"],
                call["confidence"],
                call["method"],
            )
        else:
            # no genus resolved -> best non-unclassified by rank, else unclassified
            ranked = [
                c for c in per_gene.values() if c["genus_call"] != tlca.UNCLASSIFIED
            ]
            call = (
                ranked[0]
                if ranked
                else {
                    "genus_call": tlca.UNCLASSIFIED,
                    "rank": "none",
                    "confidence": "0.000",
                    "method": "lca",
                }
            )
            genus_call, rank, conf, method = (
                call["genus_call"],
                call["rank"],
                call["confidence"],
                call["method"],
            )
        # mosaic only over main genes (exclude OTHER; ENV noisy but kept as a main gene)
        distinct = {c["genus_call"] for g, c in confident.items() if g in main_set}
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
                "completeness": f"{len(present_main) / len(main_probes):.3f}"
                if main_probes
                else "0.000",
                "canonical_order": str(canonical),
                "genus_call": genus_call,
                "rank": rank,
                "confidence": conf,
                "method": method,
                "per_gene": ";".join(
                    f"{g}:{c['genus_call']}({c['method']},{c['confidence']})"
                    for g, c in sorted(per_gene.items())
                ),
                "is_mosaic": str(len(distinct) > 1),
                "mosaic_composition": (
                    ";".join(
                        f"{g}:{c['genus_call']}" for g, c in sorted(confident.items())
                    )
                    if len(distinct) > 1
                    else ""
                ),
                "erv_class": tlca.ERV_CLASS.get(genus_call, ""),
                "probe_label_set": lc["probe_label_set"],
                "ref_version": ref_version,
            }
        )
    return records


def summarise(records: list[dict[str, str]]) -> str:
    total = len(records)
    placed = [r for r in records if r["genus_call"] != tlca.UNCLASSIFIED]
    by_rank = Counter(r["rank"] for r in placed)
    genera = Counter(r["genus_call"] for r in placed if r["rank"] == "genus")
    methods = Counter(r["method"] for r in placed if r["rank"] == "genus")
    mosaic = sum(r["is_mosaic"] == "True" for r in records)
    lines = [
        f"loci: {total}   classified: {len(placed)}   mosaic: {mosaic}",
        "rank: " + ", ".join(f"{k}={v}" for k, v in by_rank.items()),
        "genus-call method: " + ", ".join(f"{k}={v}" for k, v in methods.items()),
        "confident genera:",
    ]
    for g, n in genera.most_common():
        lines.append(f"  {g:18s} {n:6d}  [{tlca.ERV_CLASS.get(g, '')}]")
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
    "genus_call",
    "rank",
    "confidence",
    "method",
    "per_gene",
    "is_mosaic",
    "mosaic_composition",
    "erv_class",
    "probe_label_set",
    "ref_version",
]


def write_tables(records: list[dict[str, str]], parquet: Path, csv_path: Path) -> None:
    """Write the per-locus genus-call table as dual parquet + csv (always, even if empty)."""
    df = pd.DataFrame(records, columns=LOCI_COLUMNS)
    parquet.parent.mkdir(parents=True, exist_ok=True)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(parquet, index=False)
    df.to_csv(csv_path, index=False)


def write_track(records: list[dict[str, str]], gff3: Path, bed: Path) -> None:
    """Project per-locus calls to genome coordinates: IGV GFF3 + BED (colour-by-genus)."""
    gff3.parent.mkdir(parents=True, exist_ok=True)
    with gff3.open("w", encoding="utf-8") as g, bed.open("w", encoding="utf-8") as b:
        g.write("##gff-version 3\n")
        for r in records:
            attrs = (
                f"ID={r['id']};genus={r['genus_call']};rank={r['rank']};"
                f"method={r['method']};confidence={r['confidence']};"
                f"mosaic={r['is_mosaic']};erv_class={r['erv_class']};"
                f"genes={r['genes_present']}"
            )
            g.write(
                f"{r['seqname']}\tRetroSeek\tERV_locus\t{r['start']}\t{r['end']}\t"
                f"{r['confidence']}\t{r['strand']}\t.\t{attrs}\n"
            )
            # BED is 0-based half-open; name carries the genus call for IGV colour-by-name.
            b.write(
                f"{r['seqname']}\t{int(r['start']) - 1}\t{r['end']}\t"
                f"{r['id']}|{r['genus_call']}\t0\t{r['strand']}\n"
            )


def main() -> int:
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
    a = p.parse_args()

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
    )
    print(summarise(records))

    if a.out and records:  # ad-hoc single CSV (back-compat for trial docs)
        with a.out.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(records[0].keys()))
            w.writeheader()
            w.writerows(records)
        print(f"\nwrote -> {a.out}", file=sys.stderr)
    if a.out_parquet and a.out_csv:  # production dual-table output
        write_tables(records, a.out_parquet, a.out_csv)
        print(f"wrote tables -> {a.out_parquet}, {a.out_csv}", file=sys.stderr)
    if a.out_gff3 and a.out_bed:  # production IGV track
        write_track(records, a.out_gff3, a.out_bed)
        print(f"wrote track -> {a.out_gff3}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
