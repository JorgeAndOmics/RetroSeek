"""Annotate native solo-LTR calls with RetroSeek's per-locus taxonomy.

What a solo LTR is here
-----------------------
When a provirus's two flanking LTRs recombine homologously, the internal region
is excised and a single LTR is left at the integration site. Each solo LTR marks
one ancestral integration whose provirus is gone, so they extend the ERV
inventory into events the probe-based search cannot see - there is no retroviral
protein left at the locus to find.

Where they come from
--------------------
``solo_finder.py`` (ADR-017), which blasts the LTR arms of ERV-bearing elements
against the genome and subtracts the hits that are not solos. Its output is six
tab-separated columns::

    chrom    start    end    chrom:start..end    library_id    coverage

This module was written for LTR_retriever's ``solo_finder.pl``, which emitted the
same six columns, and is ported unchanged apart from this note: the six-column
table is the seam between detection and annotation, so replacing the detector
behind it required no change here. That is also why its tests carried over intact.

How a solo inherits a taxon
---------------------------
The detector names each solo after the genomic span of the element whose LTR arm
caught it, in the shape ``{chr}:{start}..{end}#LTR/{element}``. The library ID is
therefore a **coordinate**, not an opaque ``family1``, which makes the mapping
back to RetroSeek a plain interval join:

1. **Primary path - library element.** Parse the solo's ``library_id`` into a
   span and overlap it against the classified LTR-flanked loci
   (``{genome}.loci.csv``). The solo inherits ``taxon_call`` / ``rank`` /
   ``segment`` / ``erv_class`` from the locus with the largest overlap.
   ``label_source=library``. This follows sequence homology: the solo's DNA
   matched *that* element's library entry, wherever the two sit on the
   chromosome.
2. **Fallback - nearest locus.** If the library ID carries no coordinates, or
   its span overlaps no classified locus, inherit from the nearest classified
   locus on the same chromosome within ``max_distance`` bp.
   ``label_source=nearest_locus``. Proximity is a weaker signal than homology,
   so it is recorded distinctly and downstream analyses can filter on it.
3. Otherwise ``label_source=none`` and the taxonomy fields stay empty. A solo is
   never dropped for being unlabelled - the count is a real observation even
   when its lineage is not resolvable.

Overlap rather than exact coordinate equality is deliberate: a RetroSeek locus
sits *inside* its LTR element rather than sharing its edges, since the locus is
where the probe hit and the element is the whole LTR-to-LTR span.

Only ``source=ltr-flanked`` loci act as donors. An orphan has no LTR element by
definition, so it cannot have contributed a bait arm. (``ltr-flanked`` here is the
catalog's source *value*, which ADR-016 left unchanged when it renamed the track
directory from ``valid`` to ``element_hits``.)

Outputs
-------
``{genome}.gff3``
    The solo-LTR track, one feature per solo, taxonomy in the attributes.
``{genome}.solo_ltr.csv`` / ``.parquet``
    Per-solo table in ``catalog.csv``'s column vocabulary, so solos can join the
    catalog as the ``solo-ltr`` tier beside ``ltr-flanked`` and ``orphan``.
``{genome}.csv`` / ``.parquet`` (ratio)
    Solo/intact counts and ratio per group, where the denominator is the count
    of LTR-flanked loci in that group. Grouping follows ADR-012's vocabulary:
    ``segment`` (default), ``taxon_call`` or ``none``.

Interpreting the ratio
----------------------
A high solo/intact ratio means many of a lineage's integrations have had time to
recombine away, so the ratio orders lineages by relative age. It is a crude
proxy: solo LTRs degrade faster than intact elements, integration preferences
differ between lineages, and some proviruses are lost by deletion rather than
recombination. The denominator counts proviruses with retroviral *gene*
evidence, so it is also bounded by probe coverage.
"""

from __future__ import annotations

import argparse
import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from log import OK, job_logging, run_main
from tabular import tab_rows

logger = logging.getLogger(__name__)

VALID_GROUP_BY = ("segment", "taxon_call", "none")

# >{chr}:{start}..{end}[_REGION][#LTR/{family}] - annotate_lib.pl. RepeatMasker
# reports the name without the #class suffix, so that part is optional.
_LIBRARY_ID_RE = re.compile(r"^(?P<chrom>.+?):(?P<start>\d+)\.\.(?P<end>\d+)")


@dataclass
class ClassifiedLocus:
    """One classified LTR-flanked locus - a potential taxonomy donor."""

    id: str
    seqname: str
    start: int
    end: int
    taxon_call: str
    rank: str
    segment: str
    segment_rank: str
    erv_class: str
    confidence: str

    def overlap_with(self, seqname: str, start: int, end: int) -> int:
        """Return the number of shared bases, 0 when disjoint or on another chrom."""
        if self.seqname != seqname:
            return 0
        return max(0, min(self.end, end) - max(self.start, start) + 1)

    def distance_to(self, seqname: str, start: int, end: int) -> int | None:
        """Minimum bp gap on the same chromosome; 0 on overlap, None if different."""
        if self.seqname != seqname:
            return None
        if self.end < start:
            return start - self.end
        if self.start > end:
            return self.start - end
        return 0


@dataclass
class SoloLTR:
    """One solo LTR from ``solo_finder.pl``, with inherited taxonomy."""

    chrom: str
    start: int
    end: int
    library_id: str
    coverage: float
    taxon_call: str = ""
    rank: str = ""
    segment: str = ""
    segment_rank: str = ""
    erv_class: str = ""
    confidence: str = ""
    source_loci: list[str] = field(default_factory=list)
    label_source: str = "none"  # "library" | "nearest_locus" | "none"


# ---------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------
def parse_library_locus(library_id: str) -> tuple[str, int, int] | None:
    """Return ``(chrom, start, end)`` encoded in an LTR library sequence name.

    **Coordinates are normalised so start <= end.** LTR_retriever writes
    minus-strand elements in descending order (``CM040301.1:40850808..40843686``),
    and on the Desmodus pilot that was 47% of all solos - every one of them
    silently failed the overlap join and lost its taxon, because a reversed
    interval makes the shared-base arithmetic return zero rather than error.

    Returns ``None`` when the name carries no coordinates - an older
    LTR_retriever naming scheme, or a RepeatMasker rename. That is not an error;
    the caller falls back to the nearest classified locus.
    """
    match = _LIBRARY_ID_RE.match(library_id)
    if match is None:
        return None
    start, end = int(match.group("start")), int(match.group("end"))
    return match.group("chrom"), min(start, end), max(start, end)


def _solo_from_fields(fields: list[str]) -> SoloLTR | None:
    """One solo from a six-column row, or None when a number does not parse."""
    try:
        start, end = int(fields[1]), int(fields[2])
        coverage = float(fields[5])
    except ValueError:
        return None
    return SoloLTR(
        chrom=fields[0], start=start, end=end, library_id=fields[4], coverage=coverage
    )


def parse_solo_list(path: Path) -> list[SoloLTR]:
    """Read ``solo_finder.pl`` output. A missing or empty file yields no solos.

    Zero solos is a legitimate result for a genome, so absence is not an error
    here - but ``run_ltr_retriever.py`` does fail loudly when the RepeatMasker
    table that feeds solo_finder is missing, which is a different claim. Rows
    whose coordinates or coverage do not parse are skipped with a warning: the
    detector never writes one, so they mean its column layout changed.
    """
    if not path.exists():
        return []
    with path.open() as handle:
        parsed = [_solo_from_fields(fields) for fields in tab_rows(handle, 6)]
    solos = [solo for solo in parsed if solo is not None]
    if len(solos) < len(parsed):
        logger.warning(
            "%s: %d of %d solo rows have unreadable coordinates or coverage and "
            "were skipped; check the solo list's column layout",
            path.name,
            len(parsed) - len(solos),
            len(parsed),
        )
    return solos


def parse_loci_csv(path: Path) -> list[ClassifiedLocus]:
    """Read the classified LTR-flanked loci that can donate a taxon.

    Orphan rows are skipped: an orphan has no LTR element, so it cannot have
    seeded a sequence in LTR_retriever's LTR library.
    """
    if not path.exists():
        raise FileNotFoundError(f"classified loci CSV not found: {path}")
    frame = pd.read_csv(path, dtype=str).fillna("")
    if "source" in frame.columns:
        frame = frame[frame["source"] == "ltr-flanked"]
    loci: list[ClassifiedLocus] = []
    for row in frame.to_dict("records"):
        try:
            start, end = int(row["start"]), int(row["end"])
        except (KeyError, ValueError):
            continue
        loci.append(
            ClassifiedLocus(
                id=str(row.get("id", "")),
                seqname=str(row.get("seqname", "")),
                start=start,
                end=end,
                taxon_call=str(row.get("taxon_call", "")),
                rank=str(row.get("rank", "")),
                segment=str(row.get("segment", "")),
                segment_rank=str(row.get("segment_rank", "")),
                erv_class=str(row.get("erv_class", "")),
                confidence=str(row.get("confidence", "")),
            )
        )
    return loci


# ---------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------
def _inherit(solo: SoloLTR, donor: ClassifiedLocus, label_source: str) -> None:
    solo.taxon_call = donor.taxon_call
    solo.rank = donor.rank
    solo.segment = donor.segment
    solo.segment_rank = donor.segment_rank
    solo.erv_class = donor.erv_class
    solo.confidence = donor.confidence
    solo.label_source = label_source


def _library_donors(
    solo: SoloLTR, by_chrom: dict[str, list[ClassifiedLocus]]
) -> list[ClassifiedLocus]:
    """Loci overlapping the element named by the solo's library id.

    Widest overlap first, locus id breaking ties. Empty when the id carries no
    coordinates or its span overlaps no classified locus.
    """
    span = parse_library_locus(solo.library_id)
    if span is None:
        return []
    chrom, start, end = span
    widths = (
        (locus.overlap_with(chrom, start, end), locus)
        for locus in by_chrom.get(chrom, [])
    )
    hits = sorted(
        ((width, locus) for width, locus in widths if width > 0),
        key=lambda pair: (-pair[0], pair[1].id),
    )
    return [locus for _, locus in hits]


def _nearest_donor(
    solo: SoloLTR, by_chrom: dict[str, list[ClassifiedLocus]], max_distance: int
) -> ClassifiedLocus | None:
    """The closest classified locus on the solo's chromosome within reach.

    On a tie in distance the locus listed first wins (strict ``<`` below).
    """
    nearest: ClassifiedLocus | None = None
    nearest_distance = max_distance + 1
    for locus in by_chrom.get(solo.chrom, []):
        distance = locus.distance_to(solo.chrom, solo.start, solo.end)
        if distance is not None and distance < nearest_distance:
            nearest, nearest_distance = locus, distance
    return nearest


def annotate_solos(
    solos: list[SoloLTR],
    loci: list[ClassifiedLocus],
    max_distance: int,
) -> None:
    """Populate each solo's taxonomy in place (library path, then nearest locus)."""
    by_chrom: dict[str, list[ClassifiedLocus]] = defaultdict(list)
    for locus in loci:
        by_chrom[locus.seqname].append(locus)

    for solo in solos:
        donors = _library_donors(solo, by_chrom)
        if donors:
            _inherit(solo, donors[0], "library")
            solo.source_loci = [locus.id for locus in donors]
            continue
        nearest = _nearest_donor(solo, by_chrom, max_distance)
        if nearest is not None:
            _inherit(solo, nearest, "nearest_locus")
            solo.source_loci = [nearest.id]


# ---------------------------------------------------------------------
# Ratio
# ---------------------------------------------------------------------
def _group_value(taxon_call: str, segment: str, group_by: str) -> str:
    if group_by == "none":
        return "all"
    return (segment if group_by == "segment" else taxon_call) or "unassigned"


def compute_solo_intact_ratio(
    solos: list[SoloLTR],
    loci: list[ClassifiedLocus],
    species: str,
    group_by: str = "segment",
) -> pd.DataFrame:
    """Return per-group solo and intact counts with their ratio.

    Groups with intact loci but no solos are emitted with ``solo_count = 0``:
    a lineage whose proviruses have not recombined away is a finding, and
    dropping the row would silently turn it into a missing lineage.
    """
    if group_by not in VALID_GROUP_BY:
        raise ValueError(
            f"Unknown group_by {group_by!r}; expected one of {VALID_GROUP_BY}"
        )

    intact: dict[str, int] = defaultdict(int)
    for locus in loci:
        intact[_group_value(locus.taxon_call, locus.segment, group_by)] += 1
    solo: dict[str, int] = defaultdict(int)
    for entry in solos:
        solo[_group_value(entry.taxon_call, entry.segment, group_by)] += 1

    rows = []
    for group in sorted(set(intact) | set(solo)):
        solo_count, intact_count = solo.get(group, 0), intact.get(group, 0)
        total = solo_count + intact_count
        rows.append(
            {
                "species": species,
                "group": group,
                "group_by": group_by,
                "intact_count": intact_count,
                "solo_count": solo_count,
                "total_integrations": total,
                "solo_to_intact_ratio": (
                    solo_count / intact_count if intact_count else float("nan")
                ),
                "solo_fraction": (solo_count / total) if total else 0.0,
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------
def solo_id(index: int) -> str:
    """Stable per-solo identifier, shared by the GFF3 track and the table.

    The two outputs describe the same records, so they must agree: without a
    common key there is no way to look up a track feature's row in the catalog.
    Unique within a genome, like the classifier's ``L``-prefixed locus ids.
    """
    return f"S{index}"


def write_solo_ltr_gff3(solos: list[SoloLTR], output_path: Path, genome: str) -> None:
    """Write the solo-LTR track. An empty set still emits a valid header."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        handle.write("##gff-version 3\n")
        handle.write(f"##source-version RetroSeek solo_ltr_integrator {genome}\n")
        for index, solo in enumerate(solos):
            attributes = ";".join(
                (
                    f"ID={solo_id(index)}",
                    f"taxon_call={solo.taxon_call or 'unassigned'}",
                    f"segment={solo.segment or 'unassigned'}",
                    f"erv_class={solo.erv_class or 'NA'}",
                    f"library_id={solo.library_id}",
                    f"coverage={solo.coverage:.3f}",
                    f"source_loci={','.join(solo.source_loci) or 'none'}",
                    f"label_source={solo.label_source}",
                )
            )
            handle.write(
                "\t".join(
                    (
                        solo.chrom,
                        "LTR_retriever",
                        "solo_LTR",
                        str(solo.start),
                        str(solo.end),
                        ".",
                        ".",
                        ".",
                        attributes,
                    )
                )
                + "\n"
            )


def solo_table(solos: list[SoloLTR], species: str) -> pd.DataFrame:
    """Return the per-solo table in ``catalog.csv``'s column vocabulary.

    A solo LTR has no internal region by definition, so the gene-content columns
    are constants: ``structure_class=solo_ltr``, no main genes, zero
    completeness. ``source=solo-ltr`` is the third catalog tier.
    """
    return pd.DataFrame(
        [
            {
                "species": species,
                "source": "solo-ltr",
                "seqname": solo.chrom,
                "start": solo.start,
                "end": solo.end,
                "strand": ".",
                "taxon_call": solo.taxon_call,
                "rank": solo.rank,
                "segment": solo.segment,
                "segment_rank": solo.segment_rank,
                "resolved": bool(solo.taxon_call),
                "confidence": solo.confidence,
                "erv_class": solo.erv_class,
                "structure_class": "solo_ltr",
                "domain_tier": "non_domain",
                "completeness": 0.0,
                "n_main_genes": 0,
                "genes_present": "",
                "is_mosaic": False,
                "coverage": solo.coverage,
                "library_id": solo.library_id,
                "source_loci": ",".join(solo.source_loci),
                "label_source": solo.label_source,
                "id": solo_id(index),
            }
            for index, solo in enumerate(solos)
        ],
        columns=[
            "species",
            "source",
            "seqname",
            "start",
            "end",
            "strand",
            "taxon_call",
            "rank",
            "segment",
            "segment_rank",
            "resolved",
            "confidence",
            "erv_class",
            "structure_class",
            "domain_tier",
            "completeness",
            "n_main_genes",
            "genes_present",
            "is_mosaic",
            "coverage",
            "library_id",
            "source_loci",
            "label_source",
            "id",
        ],
    )


def write_solo_table(
    solos: list[SoloLTR], csv_path: Path, parquet_path: Path, species: str
) -> None:
    """Write the per-solo table as both user-facing CSV and pipeline parquet."""
    frame = solo_table(solos, species)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    parquet_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(csv_path, index=False)
    frame.to_parquet(parquet_path, index=False)


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------
def main(argv: list[str] | None = None) -> None:
    """Entry point."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--solo-list",
        type=Path,
        required=True,
        help="solo_finder.py output: the six-column solo list.",
    )
    parser.add_argument(
        "--loci-csv",
        type=Path,
        required=True,
        help="Classified loci ({genome}.loci.csv) - taxonomy donors.",
    )
    parser.add_argument("--genome", required=True, help="Genome stem, for messages.")
    parser.add_argument(
        "--species", required=True, help="Display species name, keyed to the catalog."
    )
    parser.add_argument("--output-gff3", type=Path, required=True)
    parser.add_argument("--output-solo-csv", type=Path, required=True)
    parser.add_argument("--output-solo-parquet", type=Path, required=True)
    parser.add_argument("--output-ratio-csv", type=Path, required=True)
    parser.add_argument("--output-ratio-parquet", type=Path, required=True)
    parser.add_argument("--group-by", choices=VALID_GROUP_BY, default="segment")
    parser.add_argument("--nearest-locus-max-distance", type=int, default=10000)
    parser.add_argument("--log", type=Path, help="job log (the Snakemake log: path)")
    args = parser.parse_args(argv)
    job_logging(args.log, "solo_annotator")

    loci = parse_loci_csv(args.loci_csv)
    solos = parse_solo_list(args.solo_list)
    annotate_solos(solos, loci, max_distance=args.nearest_locus_max_distance)

    by_source = dict.fromkeys(("library", "nearest_locus", "none"), 0)
    for solo in solos:
        by_source[solo.label_source] += 1
    if loci and not solos:
        logger.warning(
            "no solo LTR despite %d LTR-flanked loci; in mammals solos normally "
            "outnumber intact proviruses. Check this genome's solo funnel "
            "(tables/solo_ltr/<genome>.funnel.csv) for the step that lost them",
            len(loci),
        )

    write_solo_ltr_gff3(solos, args.output_gff3, genome=args.genome)
    write_solo_table(
        solos, args.output_solo_csv, args.output_solo_parquet, species=args.species
    )
    ratio = compute_solo_intact_ratio(
        solos, loci, species=args.species, group_by=args.group_by
    )
    args.output_ratio_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_ratio_parquet.parent.mkdir(parents=True, exist_ok=True)
    ratio.to_csv(args.output_ratio_csv, index=False)
    ratio.to_parquet(args.output_ratio_parquet, index=False)
    logger.log(
        OK,
        "%s solo LTRs annotated (taxon from library %s, nearest locus %s, "
        "unassigned %s)",
        f"{len(solos):,}",
        f"{by_source['library']:,}",
        f"{by_source['nearest_locus']:,}",
        f"{by_source['none']:,}",
    )


if __name__ == "__main__":
    run_main(main)
