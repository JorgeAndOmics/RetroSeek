"""Call solo LTRs from the bait blastn hits, by subtraction (ADR-017).

A solo LTR is what a provirus leaves behind when its two LTRs recombine
homologously and excise everything between them. There is no sequence feature that
says "solo": the evidence is genomic context, so the detector finds every copy of a
known-retroviral LTR and then removes the copies that are demonstrably something
else.

Every LTR in a genome began as one of a pair flanking a provirus, which leaves
exactly three fates and makes the subtraction complete:

    overlaps an intact LTRharvest element  -> one flank of a catalogued element
    close to an orphan locus               -> monoLTR beside surviving coding
                                              sequence, so the provirus is damaged
                                              rather than excised
    neither                                -> solo LTR

The middle class is the one LTR_retriever structurally cannot represent, because it
subtracts only its own intact elements and has no concept of an orphan. RetroSeek's
orphan track is precisely a map of retroviral coding sequence that LTRharvest
missed, which is what makes this detector possible.

Acceptance criteria come from Ou and Jiang 2018 as published (coverage window,
absolute alignment minimum) plus one of ours (bait and hit both near full length at
high identity). Ours is not optional: their coverage rule is a fraction of a
full-length family *consensus*, while our bait is individual arms as short as
100 bp, and without the length floor the method returns a 492:1 solo/intact ratio,
an order of magnitude above anything in the literature.

Outputs:
    the solo list, a six-column TSV that `solo_annotator.py` reads;
    the full candidate table, all three classes with their evidence, because
    evidence is recorded rather than gated (the ADR-015 principle);
    the funnel counts, which say where candidates died.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from collections import Counter, defaultdict
from collections.abc import Iterator
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import NamedTuple, TextIO

import pandas as pd
from solo_intervals import IntervalIndex

# The three fates. Module constants rather than string literals at call sites, so a
# typo is an ImportError instead of a silently empty class.
INTACT_FLANK = "intact_flank"
MONO_AT_ORPHAN = "mono_ltr_at_orphan"
SOLO = "solo"

ELEMENT_FEATURE = "LTR_retrotransposon"


class Hit(NamedTuple):
    """One blastn HSP, in the column order the search is asked to emit."""

    bait: str
    seqname: str
    identity: float
    alignment_length: int
    bait_length: int
    subject_start: int
    subject_end: int

    @property
    def coverage(self) -> float:
        """Alignment length as a fraction of the bait arm, Ou and Jiang's measure."""
        return self.alignment_length / self.bait_length if self.bait_length else 0.0

    @property
    def span(self) -> tuple[int, int]:
        """Subject coordinates, low first: blastn reports them reversed on minus."""
        if self.subject_start <= self.subject_end:
            return self.subject_start, self.subject_end
        return self.subject_end, self.subject_start


@dataclass(frozen=True)
class Thresholds:
    """Acceptance criteria, all supplied by the caller from config.yaml.

    No defaults: a default here would be a second source of truth that silently
    drifts from config.yaml the first time one of these changes.
    """

    min_identity: float
    min_coverage: float
    max_coverage: float
    min_alignment_length: int
    min_hit_length: int
    merge_gap: int
    orphan_pad: int


@dataclass
class Candidate:
    """A merged run of accepted hits at one genomic location."""

    seqname: str
    start: int
    end: int
    best_identity: float
    n_hits: int
    bait: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class ClassifiedCandidate:
    """A candidate with its fate and the evidence behind it."""

    candidate: Candidate
    fate: str
    orphan_distance: int | None


@dataclass
class Result:
    candidates: list[ClassifiedCandidate]
    funnel: Counter[str]


def rejection_reason(hit: Hit, thresholds: Thresholds) -> str | None:
    """Why this hit is not acceptable, or None if it is.

    Returning the reason rather than a bool is what lets the funnel plot show where
    candidates die. Order matters only for the tally, not the verdict: the absolute
    alignment floor is checked first because it is the cheapest and the most basic.
    Every comparison is inclusive, so a value exactly at a threshold passes.
    """
    if hit.alignment_length < thresholds.min_alignment_length:
        return "alignment_length"
    if not thresholds.min_coverage <= hit.coverage <= thresholds.max_coverage:
        return "coverage"
    if hit.alignment_length < thresholds.min_hit_length:
        return "hit_length"
    if hit.identity < thresholds.min_identity:
        return "identity"
    return None


def parse_hits(handle: TextIO) -> Iterator[Hit]:
    """Read the tabular blastn output, skipping rows too short to interpret."""
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 7:
            continue
        try:
            yield Hit(
                bait=fields[0],
                seqname=fields[1],
                identity=float(fields[2]),
                alignment_length=int(fields[3]),
                bait_length=int(fields[4]),
                subject_start=int(fields[5]),
                subject_end=int(fields[6]),
            )
        except ValueError:
            continue


def merge_candidates(accepted: list[Hit], gap: int) -> list[Candidate]:
    """Merge accepted hits into candidate loci, per sequence.

    Two hits that overlap are two alignments of the same genomic LTR, usually from
    the two arms of the same element, so they must count once. The representative
    bait is the best-matching one, because that is the element the solo will inherit
    its taxonomy from.
    """
    by_sequence: dict[str, list[Hit]] = defaultdict(list)
    for hit in accepted:
        by_sequence[hit.seqname].append(hit)

    candidates: list[Candidate] = []
    for seqname in sorted(by_sequence):
        hits = sorted(by_sequence[seqname], key=lambda h: h.span)
        current: Candidate | None = None
        for hit in hits:
            start, end = hit.span
            if current is not None and start <= current.end + gap:
                current.end = max(current.end, end)
                current.n_hits += 1
                if hit.identity > current.best_identity:
                    current.best_identity = hit.identity
                    current.bait = hit.bait
                continue
            if current is not None:
                candidates.append(current)
            current = Candidate(seqname, start, end, hit.identity, 1, hit.bait)
        if current is not None:
            candidates.append(current)
    return candidates


def classify(
    candidate: Candidate,
    elements: dict[str, IntervalIndex],
    orphans: dict[str, IntervalIndex],
    orphan_pad: int,
) -> str:
    """Decide a candidate's fate. Element overlap is checked first and wins.

    The priority is not arbitrary. An intact element with an orphan nearby is a
    catalogued element; calling it a monoLTR would remove it from the solo/intact
    ratio's denominator and inflate the ratio.
    """
    element_index = elements.get(candidate.seqname)
    if element_index is not None and element_index.overlaps(
        candidate.start, candidate.end
    ):
        return INTACT_FLANK
    orphan_index = orphans.get(candidate.seqname)
    if orphan_index is not None and orphan_index.overlaps(
        candidate.start, candidate.end, pad=orphan_pad
    ):
        return MONO_AT_ORPHAN
    return SOLO


def read_intervals(gff3: Path, feature: str | None = None) -> dict[str, IntervalIndex]:
    """Build a per-sequence interval index from a GFF3.

    `feature` None means every feature line counts, which is what the orphan track
    needs: its locus lines carry the probe type, not a fixed feature name.
    """
    raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with gff3.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            if feature is not None and fields[2] != feature:
                continue
            try:
                raw[fields[0]].append((int(fields[3]), int(fields[4])))
            except ValueError:
                continue
    return {seqname: IntervalIndex(spans) for seqname, spans in raw.items()}


def element_spans(gff3: Path) -> dict[str, tuple[str, int, int]]:
    """Element ID to its genomic span, for naming the seeding element of a solo.

    The integrator maps a solo back to a classified locus by parsing coordinates out
    of the library id, so the span has to travel with the call.
    """
    spans: dict[str, tuple[str, int, int]] = {}
    with gff3.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != ELEMENT_FEATURE:
                continue
            for attribute in fields[8].rstrip().split(";"):
                key, _, value = attribute.partition("=")
                if key.strip() == "ID" and value:
                    spans[value.strip()] = (fields[0], int(fields[3]), int(fields[4]))
    return spans


def _open_hits(path: Path) -> TextIO:
    """Hit tables are gzipped: 421 MB raw for one bat genome."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def run(
    hits_path: Path,
    elements_gff3: Path,
    orphans_gff3: Path,
    thresholds: Thresholds,
) -> Result:
    """Stream the hits, merge what survives, and classify each candidate."""
    funnel: Counter[str] = Counter()
    accepted: list[Hit] = []
    with _open_hits(hits_path) as handle:
        for hit in parse_hits(handle):
            funnel["raw_hits"] += 1
            reason = rejection_reason(hit, thresholds)
            if reason:
                funnel[f"rejected_{reason}"] += 1
                continue
            funnel["accepted_hits"] += 1
            accepted.append(hit)

    candidates = merge_candidates(accepted, thresholds.merge_gap)
    funnel["merged_candidates"] = len(candidates)

    elements = read_intervals(elements_gff3, ELEMENT_FEATURE)
    orphans = read_intervals(orphans_gff3)

    classified: list[ClassifiedCandidate] = []
    for candidate in candidates:
        fate = classify(candidate, elements, orphans, thresholds.orphan_pad)
        funnel[fate] += 1
        orphan_index = orphans.get(candidate.seqname)
        distance = (
            orphan_index.nearest_distance(candidate.start, candidate.end)
            if orphan_index is not None
            else None
        )
        classified.append(ClassifiedCandidate(candidate, fate, distance))
    return Result(classified, funnel)


def write_solo_list(
    classified: list[ClassifiedCandidate],
    spans: dict[str, tuple[str, int, int]],
    path: Path,
) -> int:
    """Write the six-column TSV the integrator reads, solos only.

    Columns are `chrom start end span library_id coverage`. The library id carries
    the seeding element's span in LTR_retriever's own `{chrom}:{start}..{end}#LTR/
    {family}` shape, which is what lets the integrator resolve a solo to a
    classified locus by coordinate without knowing anything about this detector.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with path.open("w") as handle:
        for item in classified:
            if item.fate != SOLO:
                continue
            candidate = item.candidate
            parent = (
                candidate.bait.split("|")[1]
                if "|" in candidate.bait
                else candidate.bait
            )
            span = spans.get(parent)
            library_id = (
                f"{span[0]}:{span[1]}..{span[2]}#LTR/{parent}"
                if span
                else f"unknown#LTR/{parent}"
            )
            handle.write(
                f"{candidate.seqname}\t{candidate.start}\t{candidate.end}\t"
                f"{candidate.seqname}:{candidate.start}..{candidate.end}\t"
                f"{library_id}\t{candidate.best_identity / 100:.4f}\n"
            )
            written += 1
    return written


CANDIDATE_COLUMNS = (
    "seqname",
    "start",
    "end",
    "length",
    "fate",
    "best_identity",
    "n_hits",
    "bait",
    "parent",
    "orphan_distance",
)


def write_candidates(
    classified: list[ClassifiedCandidate], csv_path: Path, parquet_path: Path
) -> None:
    """Write every candidate with its class and evidence, in both formats.

    All three fates are kept, not just solos: the intact and monoLTR sets are what
    the funnel, the identity-by-class panel and the locus-level comparison with
    LTR_retriever are drawn from.
    """
    rows = []
    for item in classified:
        candidate = item.candidate
        parent = (
            candidate.bait.split("|")[1] if "|" in candidate.bait else candidate.bait
        )
        rows.append(
            {
                "seqname": candidate.seqname,
                "start": candidate.start,
                "end": candidate.end,
                "length": candidate.length,
                "fate": item.fate,
                "best_identity": candidate.best_identity,
                "n_hits": candidate.n_hits,
                "bait": candidate.bait,
                "parent": parent,
                "orphan_distance": item.orphan_distance,
            }
        )

    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(CANDIDATE_COLUMNS))
        writer.writeheader()
        writer.writerows(rows)

    parquet_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=list(CANDIDATE_COLUMNS)).to_parquet(
        parquet_path, index=False
    )


# The funnel's rows, in the order the waterfall plot reads them. Declared here so
# the plot and the table cannot disagree about what the stages are.
FUNNEL_STAGES = (
    "raw_hits",
    "rejected_alignment_length",
    "rejected_coverage",
    "rejected_hit_length",
    "rejected_identity",
    "accepted_hits",
    "merged_candidates",
    INTACT_FLANK,
    MONO_AT_ORPHAN,
    SOLO,
)


def write_funnel(funnel: Counter[str], intact_loci: int, path: Path) -> None:
    """Write the funnel counts, plus the solo/intact ratio they imply.

    The denominator is catalogued LTR-flanked ERV loci, not bait arms and not bait
    elements: the ratio compares integrations that survive as solos against
    integrations that survive intact.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["stage", "count"])
        for stage in FUNNEL_STAGES:
            writer.writerow([stage, funnel.get(stage, 0)])
        writer.writerow(["intact_loci", intact_loci])
        ratio = funnel.get(SOLO, 0) / intact_loci if intact_loci else ""
        writer.writerow(["solo_to_intact_ratio", f"{ratio:.4f}" if ratio != "" else ""])


def write_manifest(
    inputs: dict[str, Path], thresholds: Thresholds, funnel: Counter[str], path: Path
) -> None:
    """Record what ran against what, so a result can be traced to its inputs.

    Same vocabulary as the hotspot stage's manifest (generator, timestamp, input
    md5s, resolved options) so the two can be read side by side. Counts live here
    too because they are cheap and they make a stale manifest obvious.
    """

    def md5(file: Path) -> str:
        # A fingerprint for provenance, not a security boundary.
        digest = hashlib.md5()
        with file.open("rb") as handle:
            for block in iter(lambda: handle.read(1 << 20), b""):
                digest.update(block)
        return digest.hexdigest()

    lines = [
        "generator: solo_ltr/solo_finder.py",
        f"timestamp: {datetime.now(timezone.utc).isoformat()}",
        "inputs:",
    ]
    for name, file in inputs.items():
        lines.append(f"  {name}:")
        lines.append(f"    path: {file}")
        lines.append(f"    md5: {md5(file) if file.exists() else 'missing'}")
    lines.append("options:")
    for field, value in vars(thresholds).items():
        lines.append(f"  {field}: {value}")
    lines.append("counts:")
    lines.extend(f"  {stage}: {funnel.get(stage, 0)}" for stage in FUNNEL_STAGES)

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def count_intact_loci(loci_csv: Path) -> int:
    """Catalogued LTR-flanked ERV loci: the solo/intact ratio's denominator."""
    with loci_csv.open(newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--hits", type=Path, required=True, help="gzipped blastn tabular output"
    )
    parser.add_argument("--elements-gff3", type=Path, required=True)
    parser.add_argument("--orphans-gff3", type=Path, required=True)
    parser.add_argument("--loci-csv", type=Path, required=True)
    parser.add_argument("--out-solo-list", type=Path, required=True)
    parser.add_argument("--out-candidates-csv", type=Path, required=True)
    parser.add_argument("--out-candidates-parquet", type=Path, required=True)
    parser.add_argument("--out-funnel-csv", type=Path, required=True)
    parser.add_argument("--out-manifest", type=Path, required=True)
    # Thresholds, every one from config.yaml solo_ltr.*, none with a default.
    parser.add_argument("--min-identity", type=float, required=True)
    parser.add_argument("--min-coverage", type=float, required=True)
    parser.add_argument("--max-coverage", type=float, required=True)
    parser.add_argument("--min-alignment-length", type=int, required=True)
    parser.add_argument("--min-hit-length", type=int, required=True)
    parser.add_argument("--merge-gap", type=int, required=True)
    parser.add_argument("--orphan-pad", type=int, required=True)
    args = parser.parse_args(argv)

    thresholds = Thresholds(
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        max_coverage=args.max_coverage,
        min_alignment_length=args.min_alignment_length,
        min_hit_length=args.min_hit_length,
        merge_gap=args.merge_gap,
        orphan_pad=args.orphan_pad,
    )
    result = run(args.hits, args.elements_gff3, args.orphans_gff3, thresholds)
    spans = element_spans(args.elements_gff3)
    intact_loci = count_intact_loci(args.loci_csv)

    solos = write_solo_list(result.candidates, spans, args.out_solo_list)
    write_candidates(
        result.candidates, args.out_candidates_csv, args.out_candidates_parquet
    )
    write_funnel(result.funnel, intact_loci, args.out_funnel_csv)
    write_manifest(
        {
            "hits": args.hits,
            "elements_gff3": args.elements_gff3,
            "orphans_gff3": args.orphans_gff3,
            "loci_csv": args.loci_csv,
        },
        thresholds,
        result.funnel,
        args.out_manifest,
    )

    ratio = solos / intact_loci if intact_loci else float("nan")
    print(
        f"raw hits {result.funnel.get('raw_hits', 0):,} -> "
        f"accepted {result.funnel.get('accepted_hits', 0):,} -> "
        f"candidates {result.funnel.get('merged_candidates', 0):,} -> "
        f"intact {result.funnel.get(INTACT_FLANK, 0):,}, "
        f"monoLTR {result.funnel.get(MONO_AT_ORPHAN, 0):,}, "
        f"solo {solos:,}"
    )
    print(f"solo/intact = {solos:,} / {intact_loci:,} = {ratio:.1f}:1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
