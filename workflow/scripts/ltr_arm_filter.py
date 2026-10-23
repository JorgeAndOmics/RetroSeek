# =============================================================================
# ltr_arm_filter.py: drop LTR elements whose arms are LINE, SINE or DNA (ADR-026)
# =============================================================================
# LTRharvest pairs any two similar stretches into an "LTR element": two copies of
# an L1 3' end or two Alus around an ERV probe hit make the element ERV-bearing,
# and its arms then seed the solo-LTR search, which counts every L1 or Alu copy
# in the genome as a solo LTR. Measured on the model genomes (2026-09-30): 83 to
# 98% of all solos came from families whose arms Dfam calls non-LTR.
#
# This step runs after LTRdigest and before ranges_analysis. For each element that
# overlaps a probe hit (a superset of the ERV-bearing ones: about 10 to 14% more
# elements than end up in the catalog) it labels BOTH arms with Dfam and says drop
# when either arm is a non-LTR repeat. Unlabelled arms and every LTR class (ERV or
# not) stay. The labelling is per family, as in ADR-023: the arms are clustered
# with cd-hit-est at the family identity, one representative per family is searched
# with every curated Dfam model, and each arm takes its family's class. That is a
# few hundred nhmmer queries per genome, not tens of thousands.
#
# ranges_analysis.R reads the table and removes the dropped elements before it
# sorts hits into LTR-flanked and orphan, so their ERV hits become orphans.
#
# Dfam's non-LTR calls agree with UCSC RepeatMasker for 75 of 75 arms, and no
# ERV LTR in 906 was called non-LTR (measurements, 2026-10-06).
# =============================================================================

"""Label the LTR arms of probe-hit elements with Dfam; list the elements to drop."""

from __future__ import annotations

import argparse
import logging
from bisect import bisect_left, bisect_right
from collections import defaultdict
from collections.abc import Iterable, Iterator
from pathlib import Path

import pyarrow.parquet as pq
from bait_builder import ARM_FEATURE, Arm
from Bio import SeqIO
from ltr_families import ELEMENT_FEATURE, run_cdhit
from ltr_family_dfam import best_hits, read_model_classes, run_nhmmer

from log import OK, PipelineError, job_logging, run_main
from tabular import gff3_attributes, gff3_features, write_csv

logger = logging.getLogger(__name__)

VERDICT_COLUMNS = [
    "element",
    "seqname",
    "arm1_family",
    "arm1_class",
    "arm2_family",
    "arm2_class",
    "drop",
    "reason",
]
# A class label that says nothing about the arm: no model matched, or a model
# matched that carries no class. Such an arm never drops its element.
UNLABELLED = {"", "unknown"}

Hit = tuple[str, int, int]  # (sequence, start, end), 1-based, ends inclusive


def is_non_ltr(dfam_class: str | None) -> bool:
    """True for a class that is a repeat but not an LTR: LINE, SINE, DNA, satellite...

    Every class starting with LTR (ERV1, ERVK, ERVL, Gypsy, Copia, MaLR...) is kept,
    ERV or not: the filter removes elements that are not LTR elements at all, and
    leaves the ERV question to the probes.
    """
    label = (dfam_class or "").strip().lower()
    return label not in UNLABELLED and not label.startswith("ltr")


def decide(classes: Iterable[str | None]) -> tuple[bool, str]:
    """(drop, reason) for an element from the class of each of its arms.

    The element goes when either arm is a non-LTR repeat. The two arms agree in
    1,249 of 1,250 measured elements, so "either" and "both" differ by about one
    element in a thousand; "either" is the owner's rule.
    """
    for label in classes:
        if is_non_ltr(label):
            return True, f"arm is {label}"
    return False, ""


def element_arms(gff3: Path) -> list[Arm]:
    """Both arms of every LTRdigest element, numbered 1 and 2 in genomic order.

    The number is positional, not strand-aware: only the pair matters here. An
    arm row without a Parent cannot be tied to an element and is skipped.
    """
    found: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    for fields, start, end in gff3_features(gff3):
        if fields[2] != ARM_FEATURE:
            continue
        parent = gff3_attributes(fields[8]).get("Parent", "")
        if parent:
            found[parent].append((fields[0], start, end))
    return [
        Arm(seqname, start, end, parent, str(number))
        for parent, parts in found.items()
        for number, (seqname, start, end) in enumerate(
            sorted(parts, key=_by_position), start=1
        )
    ]


def _by_position(part: tuple[str, int, int]) -> tuple[int, int]:
    return part[1], part[2]


def element_spans(gff3: Path) -> dict[str, tuple[str, int, int]]:
    """{element id: (sequence, start, end)} from the LTR_retrotransposon rows."""
    spans: dict[str, tuple[str, int, int]] = {}
    for fields, start, end in gff3_features(gff3):
        if fields[2] != ELEMENT_FEATURE:
            continue
        element = gff3_attributes(fields[8]).get("ID", "")
        if element:
            spans[element] = (fields[0], start, end)
    return spans


def _index_spans(
    spans: dict[str, tuple[str, int, int]],
) -> dict[str, tuple[list[int], int, list[tuple[int, int, str]]]]:
    """Per sequence: sorted starts, the longest element, and the (start, end, id) rows."""
    by_sequence: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for element, (seqname, start, end) in spans.items():
        by_sequence[seqname].append((start, end, element))
    index = {}
    for seqname, rows in by_sequence.items():
        rows.sort()
        longest = max(end - start for start, end, _ in rows)
        index[seqname] = ([start for start, _, _ in rows], longest, rows)
    return index


def overlapping_elements(
    spans: dict[str, tuple[str, int, int]], hits: Iterable[Hit]
) -> set[str]:
    """Elements that share at least one base with any hit, whatever the strand.

    A superset of what ranges_analysis.R calls ERV-bearing: that overlap is
    strand-aware and runs on hits already filtered by bitscore and identity. A
    superset is safe, since a verdict for an element that ends up with no hit is
    never read; a subset would leave a real ERV-bearing element unlabelled.
    """
    index = _index_spans(spans)
    marked: set[str] = set()
    for seqname, first, second in hits:
        if seqname not in index:
            continue
        low, high = min(first, second), max(first, second)
        starts, longest, rows = index[seqname]
        # An element starting before low - longest cannot reach low.
        window = rows[bisect_left(starts, low - longest) : bisect_right(starts, high)]
        marked.update(element for _, end, element in window if end >= low)
    return marked


def _by_sequence(arms: list[Arm]) -> dict[str, list[Arm]]:
    grouped: dict[str, list[Arm]] = defaultdict(list)
    for arm in arms:
        grouped[arm.seqname].append(arm)
    return grouped


def _missing_sequences(remaining: set[str], fasta: Path) -> PipelineError:
    return PipelineError(
        f"{len(remaining)} sequences with LTR arms are not in {fasta}, "
        f"e.g. {sorted(remaining)[0]}",
        hint="the LTRdigest track and the genome FASTA must come from one genome",
    )


def extract_arms(fasta: Path, arms: list[Arm]) -> dict[str, str]:
    """{arm name: sequence}, cut from the genome one sequence at a time.

    Streaming keeps memory at one chromosome, not the genome. Coordinates are
    1-based with both ends inclusive, as in GFF3.

    Raises:
        PipelineError: If an arm lies on a sequence the FASTA does not have.
    """
    wanted = _by_sequence(arms)
    remaining = set(wanted)
    cut: dict[str, str] = {}
    # Biopython ships no type hints; validator.py silences the same call.
    for record in SeqIO.parse(fasta, "fasta"):  # type: ignore[no-untyped-call]
        for arm in wanted.get(record.id, []):
            cut[arm.name] = str(record.seq[arm.start - 1 : arm.end])
        remaining.discard(record.id)
        if not remaining:
            break
    if remaining:
        raise _missing_sequences(remaining, fasta)
    return cut


def verdict_rows(
    arms: list[Arm], family_of: dict[str, str], class_of: dict[str, str]
) -> list[dict[str, object]]:
    """One row per element: each arm's family and class, and the drop verdict.

    An arm that is in no family (too short to cluster) or whose family matched no
    Dfam model is unlabelled, and unlabelled arms never drop an element.
    """
    pairs: dict[str, list[Arm]] = {}
    for arm in arms:
        pairs.setdefault(arm.parent, []).append(arm)
    return [
        _verdict_row(element, pair, family_of, class_of)
        for element, pair in pairs.items()
    ]


def _verdict_row(
    element: str, pair: list[Arm], family_of: dict[str, str], class_of: dict[str, str]
) -> dict[str, object]:
    row: dict[str, object] = {"element": element, "seqname": pair[0].seqname}
    classes = []
    for number in ("1", "2"):
        arm = next((a for a in pair if a.arm == number), None)
        family = family_of.get(arm.name, "") if arm else ""
        label = class_of.get(family, "") if family else ""
        row[f"arm{number}_family"] = family
        row[f"arm{number}_class"] = label
        classes.append(label)
    drop, reason = decide(classes)
    row.update(drop=drop, reason=reason)
    return row


def read_hits(parquet: Path) -> Iterator[Hit]:
    """(accession, start, end) of every probe hit in a segmented BLAST parquet."""
    table = pq.read_table(
        parquet, columns=["accession", "hsp_sbjct_start", "hsp_sbjct_end"]
    )
    columns = (table[name].to_pylist() for name in table.column_names)
    return zip(*columns, strict=True)


def write_fasta(sequences: dict[str, str], path: Path) -> None:
    """Write `sequences` as FASTA, creating the folder."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in sequences.items()))


def family_classes(
    sequences: dict[str, str], hmm: Path, workdir: Path, identity: float, threads: int
) -> tuple[dict[str, str], dict[str, str]]:
    """(arm -> family, family -> Dfam class): cluster the arms, search one per family.

    A family is named by its representative arm. A family whose representative has
    no hit under E 1e-5 has no entry in the class table, so its arms are unlabelled.
    """
    arms_fna = workdir / "arms.fna"
    write_fasta(sequences, arms_fna)
    clusters = run_cdhit(arms_fna, workdir / "cdhit", identity, threads)
    family_of = {member: rep for rep, members in clusters for member in members}
    reps_fna = workdir / "representatives.fna"
    write_fasta({rep: sequences[rep] for rep, _ in clusters}, reps_fna)
    best = best_hits(run_nhmmer(hmm, reps_fna, workdir / "dfam", threads))
    models = read_model_classes(hmm, {hit.model for hit in best.values()})
    class_of = {rep: models.get(hit.model, "") for rep, hit in best.items()}
    return family_of, class_of


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ltrdigest-gff3", type=Path, required=True)
    parser.add_argument("--blast-parquet", type=Path, required=True)
    parser.add_argument("--genome-fasta", type=Path, required=True)
    parser.add_argument("--dfam-hmm", type=Path, required=True)
    parser.add_argument(
        "--identity",
        type=float,
        required=True,
        help="From config solo_ltr.families.identity: the arms are grouped into "
        "families exactly as the solo-LTR stage groups them.",
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--out-csv", type=Path, required=True)
    parser.add_argument("--log", type=Path, help="job log (the Snakemake log: path)")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    """Command-line entry: write the per-element verdict table for one genome."""
    args = _parse_args(argv)
    job_logging(args.log, "ltr_arm_filter")
    marked = overlapping_elements(
        element_spans(args.ltrdigest_gff3), read_hits(args.blast_parquet)
    )
    arms = [arm for arm in element_arms(args.ltrdigest_gff3) if arm.parent in marked]
    if not arms:
        logger.warning("no LTR element overlaps a probe hit: nothing to label or drop")
        write_csv([], VERDICT_COLUMNS, args.out_csv)
        return
    sequences = extract_arms(args.genome_fasta, arms)
    family_of, class_of = family_classes(
        sequences, args.dfam_hmm, args.workdir, args.identity, args.threads
    )
    rows = verdict_rows(arms, family_of, class_of)
    write_csv(rows, VERDICT_COLUMNS, args.out_csv)
    dropped = sum(bool(row["drop"]) for row in rows)
    logger.log(
        OK,
        "%s of %s probe-hit elements dropped (an arm is a non-LTR repeat in Dfam); "
        "%s arm families labelled",
        f"{dropped:,}",
        f"{len(rows):,}",
        f"{len(class_of):,}",
    )


if __name__ == "__main__":
    run_main(main)
