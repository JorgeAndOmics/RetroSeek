# =============================================================================
# ltr_families.py: group one genome's bait arms into LTR families (ADR-023)
# =============================================================================
# A family is a group of LTR arms at least `identity` identical (0.80 by default,
# the field's 80% convention) over the whole shorter arm, found by cd-hit-est.
# The arms are the solo-LTR bait: the flanking LTRs of intact, ERV-bearing
# elements. Every element joins the family of its arms, and every solo later
# inherits the family of the arm that seeded it (solo_annotator.py), as it
# inherits a genus.
#
# Measured on the model genomes (2026-09-30): at 0.80 an element's two arms land
# in one family 97% of the time and about 94% of a family's elements share one
# genus; at 0.90 and above families shatter (arm pairs together 82% or less).
#
# Identifiers are `<code>_F001`, numbered by size (largest first), ties broken by
# the representative arm's position, so the same input always gives the same
# names. The code is one genus letter plus three species letters (Mus_musculus:
# Mmus). They are unrelated to the evidence tree's `F001` clades.
# =============================================================================

"""Group one genome's solo-LTR bait arms into LTR families with cd-hit-est."""

from __future__ import annotations

import argparse
import csv
import logging
import statistics
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from bait_builder import Arm

from external import run_tool
from log import OK, PipelineError, job_logging, run_main
from tabular import gff3_attributes, gff3_features, tab_rows, write_csv
from utils import file_md5

logger = logging.getLogger(__name__)

CDHIT = "cd-hit-est"
MIN_IDENTITY = 0.80  # cd-hit-est refuses anything lower
ELEMENT_FEATURE = "LTR_retrotransposon"

FAMILY_COLUMNS = [
    "arm",
    "seqname",
    "start",
    "end",
    "element",
    "ltr_family",
    "representative",
]
SUMMARY_COLUMNS = [
    "ltr_family",
    "n_arms",
    "n_elements",
    "representative",
    "rep_seqname",
    "rep_start",
    "rep_end",
    "majority_genus",
    "genus_purity",
    "median_arm_similarity",
    "split_elements",
]
GENUS_COLUMNS = ["ltr_family", "genus", "n_elements"]

Element = tuple[str, str]  # (seqname, element id): LTRdigest numbers ids per genome


@dataclass(frozen=True)
class Family:
    """One LTR family: its identifier, its representative arm and every member arm."""

    name: str
    representative: str
    members: list[str]


# ---------------------------------------------------------------- names
def species_code(genome: str) -> str:
    """Four letters naming a genome in family identifiers: Mus_musculus gives Mmus.

    Raises:
        PipelineError: If the genome name is not `Genus_species`.
    """
    parts = genome.split("_")
    if len(parts) < 2 or not parts[0] or not parts[1]:
        raise PipelineError(
            f"cannot make a species code from the genome name {genome!r}",
            hint="name genomes Genus_species, as the config's species: block does",
        )
    return parts[0][0].upper() + parts[1][:3].lower()


def read_bait_bed(bed: Path) -> dict[str, Arm]:
    """The bait arms by name, from the bait builder's BED (0-based starts).

    The name is the bait builder's `seqname|element|L` (or `R`); the element and
    the side are read back from its last two fields.
    """
    arms = {}
    with bed.open(encoding="utf-8") as handle:
        for fields in tab_rows(handle, 4):
            _, element, side = fields[3].rsplit("|", 2)
            start, end = int(fields[1]) + 1, int(fields[2])
            arms[fields[3]] = Arm(fields[0], start, end, element, side)
    return arms


def read_genus(loci_csv: Path) -> dict[Element, str]:
    """Each classified element's genus (the `segment` column of the loci table)."""
    with loci_csv.open(newline="", encoding="utf-8") as handle:
        return {
            (r["seqname"], r["parent"]): r["segment"] for r in csv.DictReader(handle)
        }


def read_arm_similarity(ltrdigest_gff3: Path) -> dict[Element, float]:
    """Each element's LTRdigest `ltr_similarity`: how alike its two arms are, in %.

    Two arms are identical when a provirus inserts and drift apart afterwards, so
    this is an age signal.
    """
    similarity = {}
    for fields, _, _ in gff3_features(ltrdigest_gff3):
        if fields[2] != ELEMENT_FEATURE:
            continue
        attributes = gff3_attributes(fields[8])
        if "ID" in attributes and "ltr_similarity" in attributes:
            similarity[(fields[0], attributes["ID"])] = float(
                attributes["ltr_similarity"]
            )
    return similarity


# ---------------------------------------------------------------- cd-hit-est
def word_size(identity: float) -> int:
    """cd-hit-est's word length for an identity threshold (its user guide's table).

    Raises:
        PipelineError: Below 0.80, which cd-hit-est does not accept.
    """
    if identity < MIN_IDENTITY:
        raise PipelineError(
            f"LTR family identity {identity} is below 0.8, the lowest cd-hit-est takes",
            hint="set solo_ltr.families.identity between 0.8 and 1.0",
        )
    for floor, word in ((0.95, 10), (0.90, 8), (0.88, 7), (0.85, 6)):
        if identity >= floor:
            return word
    return 5


def cdhit_command(fna: Path, prefix: Path, identity: float, threads: int) -> list[str]:
    """cd-hit-est over whole arms, both strands, each arm into its best family.

    `-G 1` counts identity over the whole shorter arm; `-r 1` compares both
    strands, because bait arms are not oriented; `-g 1` puts each arm in the
    most similar family rather than the first good enough; `-d 0` keeps full
    names; `-M 0` lifts the memory cap.
    """
    # One flag and its value per line; the formatter would split each pair.
    return [
        CDHIT,
        "-i", str(fna),
        "-o", str(prefix),
        "-c", str(identity),
        "-n", str(word_size(identity)),
        "-G", "1",
        "-r", "1",
        "-g", "1",
        "-d", "0",
        "-M", "0",
        "-T", str(threads),
    ]  # fmt: skip


def parse_clstr(clstr: Path) -> list[tuple[str, list[str]]]:
    """(representative, members) per cluster of a cd-hit `.clstr` file."""
    clusters: list[tuple[str, list[str]]] = []
    for line in clstr.read_text(encoding="utf-8").splitlines():
        if line.startswith(">Cluster"):
            clusters.append(("", []))
            continue
        name = line.split(">", 1)[1].split("...")[0]
        members = clusters[-1][1]
        members.append(name)
        if line.endswith("*"):  # cd-hit marks the representative with a star
            clusters[-1] = (name, members)
    return clusters


def run_cdhit(
    fna: Path, workdir: Path, identity: float, threads: int
) -> list[tuple[str, list[str]]]:
    """Cluster the arms in `fna` and return cd-hit's clusters."""
    workdir.mkdir(parents=True, exist_ok=True)
    prefix = workdir / "families"
    run_tool(cdhit_command(fna, prefix, identity, threads))
    return parse_clstr(Path(f"{prefix}.clstr"))


def name_families(
    clusters: list[tuple[str, list[str]]], arms: dict[str, Arm], code: str
) -> list[Family]:
    """Families named `<code>_F001...`, largest first, ties by representative position.

    Raises:
        PipelineError: If a clustered arm is missing from the bait BED.
    """
    missing = sorted({m for _, members in clusters for m in members} - arms.keys())
    if missing:
        raise PipelineError(
            f"{len(missing)} clustered arms are not in the bait BED, e.g. {missing[0]}",
            hint="rerun the bait builder; the BED and FASTA must come from one run",
        )

    def order(cluster: tuple[str, list[str]]) -> tuple[int, str, int]:
        rep = arms[cluster[0]]
        return (-len(cluster[1]), rep.seqname, rep.start)

    ranked = sorted(clusters, key=order)
    width = max(3, len(str(len(ranked))))
    return [
        Family(f"{code}_F{i:0{width}d}", rep, members)
        for i, (rep, members) in enumerate(ranked, start=1)
    ]


# ---------------------------------------------------------------- tables
def family_rows(
    families: list[Family], arms: dict[str, Arm]
) -> list[dict[str, object]]:
    """One row per arm: where it sits, its element and its family."""
    return [
        {
            "arm": name,
            "seqname": arms[name].seqname,
            "start": arms[name].start,
            "end": arms[name].end,
            "element": arms[name].parent,
            "ltr_family": family.name,
            "representative": name == family.representative,
        }
        for family in families
        for name in family.members
    ]


def _elements(family: Family, arms: dict[str, Arm]) -> list[Element]:
    """The family's elements, each once, in a fixed order."""
    return sorted({(arms[m].seqname, arms[m].parent) for m in family.members})


def _majority(counts: Counter[str]) -> tuple[str, int]:
    """The most common value; a tie goes to the name first in byte order."""
    return min(counts.items(), key=lambda item: (-item[1], item[0]))


def _summary_row(
    family: Family,
    arms: dict[str, Arm],
    genus: dict[Element, str],
    similarity: dict[Element, float],
    families_of: dict[Element, set[str]],
) -> dict[str, object]:
    """One family's row of the summary table."""
    elements = _elements(family, arms)
    split = sum(len(families_of[e]) > 1 for e in elements)
    majority, count = _majority(Counter(genus.get(e, "") for e in elements))
    ages = [similarity[e] for e in elements if e in similarity]
    rep = arms[family.representative]
    return {
        "ltr_family": family.name,
        "n_arms": len(family.members),
        "n_elements": len(elements),
        "representative": family.representative,
        "rep_seqname": rep.seqname,
        "rep_start": rep.start,
        "rep_end": rep.end,
        "majority_genus": majority,
        "genus_purity": round(count / len(elements), 4),
        # LTRdigest gives two decimals; a median of two is their mean.
        "median_arm_similarity": round(statistics.median(ages), 3) if ages else "",
        "split_elements": split,
    }


def summary_rows(
    families: list[Family],
    arms: dict[str, Arm],
    genus: dict[Element, str],
    similarity: dict[Element, float],
) -> list[dict[str, object]]:
    """One row per family: size, representative, genus mix, age and split elements.

    `split_elements` counts the family's elements whose other arm fell in another
    family: the arm pair is the method's own positive control.
    """
    families_of: dict[Element, set[str]] = {}  # element -> families of its arms
    for family in families:
        for element in _elements(family, arms):
            families_of.setdefault(element, set()).add(family.name)
    return [
        _summary_row(family, arms, genus, similarity, families_of)
        for family in families
    ]


def genus_rows(
    families: list[Family], arms: dict[str, Arm], genus: dict[Element, str]
) -> list[dict[str, object]]:
    """Family by genus: how many of each family's elements carry each genus call."""
    rows: list[dict[str, object]] = []
    for family in families:
        counts = Counter(genus.get(e, "") for e in _elements(family, arms))
        rows.extend(
            {"ltr_family": family.name, "genus": g, "n_elements": n}
            for g, n in sorted(counts.items())
        )
    return rows


def write_manifest(
    inputs: dict[str, Path], identity: float, counts: dict[str, int], path: Path
) -> None:
    """Record what ran against what, as the other solo-LTR manifests do."""
    lines = [
        "generator: solo_ltr/ltr_families.py",
        f"timestamp: {datetime.now(timezone.utc).isoformat()}",
        "inputs:",
    ]
    for name, file in inputs.items():
        lines += [f"  {name}:", f"    path: {file}", f"    md5: {file_md5(file)}"]
    lines += [
        "options:",
        f"  identity: {identity}",
        f"  word_size: {word_size(identity)}",
        "  rule: identity over the whole shorter arm (cd-hit-est -G 1), both strands",
        "counts:",
    ]
    lines += [f"  {name}: {value}" for name, value in counts.items()]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------- command line
def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    """The command line: the bait, the tables it is compared with, and the outputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", required=True, help="genome name, Genus_species")
    parser.add_argument("--bait-fna", type=Path, required=True)
    parser.add_argument("--bait-bed", type=Path, required=True)
    parser.add_argument("--loci-csv", type=Path, required=True)
    parser.add_argument("--ltrdigest-gff3", type=Path, required=True)
    parser.add_argument("--identity", type=float, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--out-families", type=Path, required=True)
    parser.add_argument("--out-summary", type=Path, required=True)
    parser.add_argument("--out-genus", type=Path, required=True)
    parser.add_argument("--out-manifest", type=Path, required=True)
    parser.add_argument("--log", type=Path, help="job log (the Snakemake log: path)")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    """Command-line entry: build one genome's LTR families and write their tables."""
    args = _parse_args(argv)
    job_logging(args.log, "solo_family_builder")
    arms = read_bait_bed(args.bait_bed)
    clusters = run_cdhit(args.bait_fna, args.workdir, args.identity, args.threads)
    families = name_families(clusters, arms, species_code(args.genome))
    genus = read_genus(args.loci_csv)
    summary = summary_rows(
        families, arms, genus, read_arm_similarity(args.ltrdigest_gff3)
    )
    write_csv(family_rows(families, arms), FAMILY_COLUMNS, args.out_families)
    write_csv(summary, SUMMARY_COLUMNS, args.out_summary)
    write_csv(genus_rows(families, arms, genus), GENUS_COLUMNS, args.out_genus)
    inputs = {
        "bait_fna": args.bait_fna,
        "bait_bed": args.bait_bed,
        "loci_csv": args.loci_csv,
        "ltrdigest_gff3": args.ltrdigest_gff3,
    }
    counts = {
        "arms": len(arms),
        "families": len(families),
        "single_arm_families": sum(len(f.members) == 1 for f in families),
    }
    write_manifest(inputs, args.identity, counts, args.out_manifest)
    logger.log(OK, "%s arms in %s LTR families", f"{len(arms):,}", f"{len(families):,}")


if __name__ == "__main__":
    run_main(main)
