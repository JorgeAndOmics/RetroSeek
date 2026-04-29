"""Post-process LTR_retriever output into RetroSeek-annotated solo LTRs.

This is *Coupling B* of the LTR_retriever integration: we propagate
RetroSeek's probe labels (from ``valid_ranges.gff3``) onto the solo
LTRs that LTR_retriever discovered, using a hybrid two-tier approach
that favours LTR_retriever's own family-clustering when available and
falls back to nearest-valid-ERV location inheritance when it isn't.

Data sources
------------
``nmtf.pass.list``
    LTR_retriever's list of solo + truncated LTRs. Tab-separated. The
    first column is typically ``chrom:start..end`` or similar
    coordinate encoding; later columns include family ID, length,
    and identity. The exact format varies slightly across
    LTR_retriever releases — the parser here is defensive: it looks
    for a coordinate column and a family column by heuristic.

``LTRlib.fa``
    LTR_retriever's family consensus library. FASTA headers carry
    the family ID plus a listing of source intact-ERV IDs that seeded
    the consensus. We parse the headers to build a ``family → list of
    source-ERV-IDs`` mapping.

``valid_ranges.gff3``
    RetroSeek's domain-validated retroviral ERV track. Each row has
    genomic coordinates + ``probe`` (or ``probe_labels``) attribute
    encoding the probe family or families the ERV matched. We parse
    this into a per-chromosome interval structure so location-based
    lookup can find the valid ERV covering (or adjacent to) any
    LTR_retriever intact ERV.

Label propagation — hybrid approach
-----------------------------------
For each solo LTR in ``nmtf.pass.list``:

1. **Primary path — consensus-family mapping.** Look up the solo
   LTR's family ID in the ``family → source-ERVs`` map built from
   ``LTRlib.fa``. For each source ERV, find the corresponding valid
   range in ``valid_ranges.gff3`` by coordinate overlap. Collect the
   probe labels from those valid ranges. Label the solo LTR with the
   union.
2. **Fallback path — nearest-ERV.** If the primary path produced zero
   labels (family unresolved, source ERVs couldn't be resolved to
   valid-range entries, or the LTRlib.fa header format was
   unparseable), find the nearest valid ERV on the same chromosome
   within ``nearest_erv_max_distance`` bp. Inherit its labels.

Either way, the output GFF3 records which mechanism produced the
labels via the ``label_source`` attribute (``family`` or
``nearest_erv``) plus the list of contributing ERVs, so downstream
consumers can reaggregate as needed.

Solo / intact ratio
-------------------
In addition to the GFF3 output, this script computes a per-genome
CSV recording, per probe family (main + accessory), the solo and
intact counts and their ratio. Rows with ``label_mode='exclusive'``
count a solo/intact LTR in exactly one family (its labels were a
single-value set); rows with ``label_mode='shared'`` count
multi-labelled records in every family they claim. Both views are
emitted so downstream choice of counting rule stays flexible.

Usage
-----
::

    python solo_ltr_integrator.py \\
        --nmtf-pass-list <path> \\
        --ltr-library <path> \\
        --pass-list-gff3 <path> \\
        --valid-ranges <path> \\
        --genome <name> \\
        --output-gff3 <path> \\
        --output-ratio-csv <path> \\
        --nearest-erv-max-distance 10000
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd


# ------------------------------------------------------------------
# Data classes
# ------------------------------------------------------------------
@dataclass
class ValidRange:
    """One retroviral-confirmed ERV from valid_ranges.gff3."""

    chrom: str
    start: int  # 0-indexed inclusive
    end: int  # 0-indexed inclusive (GFF3 5th field is 1-indexed closed; we normalise)
    probes: list[str] = field(default_factory=list)
    erv_id: str = ""

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        """Closed-interval overlap test, 0-indexed."""
        return self.chrom == chrom and self.start <= end and self.end >= start

    def distance_to(self, chrom: str, start: int, end: int) -> int | None:
        """Minimum bp distance to the given interval on the same chromosome.

        Returns None if the chromosomes differ. Zero on overlap.
        """
        if self.chrom != chrom:
            return None
        if self.end < start:
            return start - self.end
        if self.start > end:
            return self.start - end
        return 0


@dataclass
class SoloLTR:
    """One solo LTR from LTR_retriever nmtf.pass.list with annotations."""

    chrom: str
    start: int  # 1-indexed (kept as emitted — we convert at GFF3 write time only)
    end: int
    strand: str = "."
    family: str | None = None
    probe_labels: list[str] = field(default_factory=list)
    contributing_ervs: list[str] = field(default_factory=list)
    label_source: str = "none"  # "family" | "nearest_erv" | "none"


# ------------------------------------------------------------------
# Parsers
# ------------------------------------------------------------------
def _parse_probes_from_gff3_attrs(attrs: str) -> list[str]:
    """Extract probe labels from a GFF3 attribute string.

    Looks for ``probe_labels=...`` (multi-label), ``probe_category=...``,
    and ``probe=...`` in that preference order. Comma-separated values
    are split into a list.
    """
    for key in ("probe_labels", "probe", "probe_category"):
        match = re.search(rf"(?:^|;){re.escape(key)}=([^;]+)", attrs)
        if match:
            raw = match.group(1).strip()
            return [p.strip() for p in raw.split(",") if p.strip()]
    return []


def parse_valid_ranges(path: Path) -> list[ValidRange]:
    """Return all ranges from valid_ranges.gff3 with their probe labels."""
    if not path.exists():
        raise FileNotFoundError(f"valid_ranges GFF3 not found: {path}")
    ranges: list[ValidRange] = []
    with path.open() as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            try:
                start = int(fields[3]) - 1
                end = int(fields[4])
            except ValueError:
                continue
            chrom = fields[0]
            attrs = fields[8]
            id_match = re.search(r"(?:^|;)ID=([^;]+)", attrs)
            erv_id = id_match.group(1) if id_match else ""
            ranges.append(
                ValidRange(
                    chrom=chrom,
                    start=start,
                    end=end,
                    probes=_parse_probes_from_gff3_attrs(attrs),
                    erv_id=erv_id,
                )
            )
    return ranges


def parse_ltr_library_headers(path: Path) -> dict[str, list[str]]:
    """Return a ``family → list of source-ERV IDs`` mapping from LTRlib.fa.

    LTR_retriever's FASTA headers look approximately like::

        >family1#LTR/unknown  members=17  source=LTR_retrotransposon5,LTR_retrotransposon12,...

    Format varies across LTR_retriever versions. The parser looks for
    a ``source=`` or ``members=`` token; if neither is present, the
    family has no traceable source ERVs and the caller's fallback
    mechanism kicks in.
    """
    if not path.exists():
        return {}
    mapping: dict[str, list[str]] = {}
    with path.open() as handle:
        for raw in handle:
            if not raw.startswith(">"):
                continue
            header = raw[1:].strip()
            # Family ID is the first whitespace-delimited token, minus
            # any ``#TE_class`` suffix.
            first = header.split()[0]
            family = first.split("#")[0]
            # Try source=<ids> then members=<ids>.
            source_match = re.search(r"(?:source|members)=([^\s]+)", header)
            source_ids: list[str] = []
            if source_match:
                source_ids = [
                    s.strip() for s in source_match.group(1).split(",") if s.strip()
                ]
            mapping[family] = source_ids
    return mapping


def _parse_coord_column(cell: str) -> tuple[str, int, int] | None:
    """Parse an LTR_retriever-style coordinate cell.

    Handles the two common formats LTR_retriever emits:

    - ``chrom:start..end``  (older releases)
    - ``chrom:start..end(+|-)``  (with strand suffix, newer)

    Returns ``(chrom, start, end)`` or ``None`` if the cell doesn't
    parse. Start / end returned as-is (LTR_retriever uses 1-indexed
    inclusive).
    """
    match = re.match(r"^([^\s:]+):(\d+)\.\.(\d+)([+\-])?\s*$", cell)
    if not match:
        return None
    return match.group(1), int(match.group(2)), int(match.group(3))


def parse_nmtf_pass_list(path: Path) -> list[SoloLTR]:
    """Return all solo / truncated LTR entries from nmtf.pass.list.

    The file is whitespace-delimited; we're defensive about column
    layout, which varies across LTR_retriever versions. We look for
    a column that parses as a coordinate tuple (``chrom:start..end``)
    and a column that looks like a family ID (``family\\d+`` or a
    bare alphanumeric family token).
    """
    if not path.exists():
        return []
    solos: list[SoloLTR] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            coord: tuple[str, int, int] | None = None
            family: str | None = None
            strand: str = "."
            for cell in cols:
                if coord is None:
                    parsed = _parse_coord_column(cell)
                    if parsed is not None:
                        coord = parsed
                        # Try to peel off a trailing strand suffix if present.
                        if cell.endswith(("+", "-")):
                            strand = cell[-1]
                        continue
                if family is None and re.match(r"^[A-Za-z][\w\-\.]+$", cell):
                    # Heuristic: a family ID-looking token AFTER the coord column.
                    # Skip a few known decoration tokens so we don't mis-identify.
                    if cell.lower() in {"solo", "truncated", "intact", "pass", "nmtf"}:
                        continue
                    family = cell
            if coord is None:
                continue
            chrom, start, end = coord
            solos.append(
                SoloLTR(
                    chrom=chrom,
                    start=start,
                    end=end,
                    strand=strand,
                    family=family,
                )
            )
    return solos


# ------------------------------------------------------------------
# Label propagation
# ------------------------------------------------------------------
def _valid_by_chrom(ranges: list[ValidRange]) -> dict[str, list[ValidRange]]:
    """Group valid ranges by chromosome for fast chromosome-scoped lookups."""
    out: dict[str, list[ValidRange]] = defaultdict(list)
    for r in ranges:
        out[r.chrom].append(r)
    for chrom in out:
        out[chrom].sort(key=lambda r: r.start)
    return dict(out)


def _resolve_source_ervs_to_valid_ranges(
    source_ids: list[str],
    valid_ranges: list[ValidRange],
) -> list[ValidRange]:
    """Match source-ERV IDs from LTR_retriever against valid_ranges.

    LTR_retriever's source IDs refer to LTRharvest candidate ERVs (e.g.
    ``LTR_retrotransposon5``) — these don't match RetroSeek's own
    probe-based IDs in ``valid_ranges.gff3``. Fortunately, LTR_retriever
    preserves the source ERV's genomic coordinates in its
    ``pass.list.gff3``. If that's parsed separately and passed here,
    we can intersect coordinates.

    For the first implementation, this helper returns an empty list
    when source IDs can't be directly matched to valid_ranges IDs —
    the caller then falls through to the nearest-ERV mechanism. Future
    enhancement: pass the parsed ``pass.list.gff3`` so we can do a
    two-step coord mapping.
    """
    # Build an ID → ValidRange map; if the ID nomenclature happens to
    # align, we resolve cleanly. Otherwise return empty and let the
    # fallback handle it.
    id_to_range = {r.erv_id: r for r in valid_ranges if r.erv_id}
    resolved: list[ValidRange] = []
    for sid in source_ids:
        if sid in id_to_range:
            resolved.append(id_to_range[sid])
    return resolved


def _nearest_valid_erv(
    solo: SoloLTR,
    valid_by_chrom: dict[str, list[ValidRange]],
    max_distance: int,
) -> ValidRange | None:
    """Return the nearest valid ERV on the same chromosome within max_distance."""
    candidates = valid_by_chrom.get(solo.chrom, [])
    if not candidates:
        return None
    # Convert solo coords to 0-indexed closed for distance calc.
    start = solo.start - 1
    end = solo.end
    best: ValidRange | None = None
    best_dist: int | None = None
    for r in candidates:
        dist = r.distance_to(solo.chrom, start, end)
        if dist is None or dist > max_distance:
            continue
        if best_dist is None or dist < best_dist:
            best_dist = dist
            best = r
    return best


def propagate_labels(
    solos: list[SoloLTR],
    valid_ranges: list[ValidRange],
    family_to_sources: dict[str, list[str]],
    max_distance: int,
) -> None:
    """In-place annotation of solo LTRs with probe labels.

    Applies the hybrid strategy: consensus-family primary path, nearest-
    ERV fallback. Each solo's ``probe_labels``, ``contributing_ervs``,
    and ``label_source`` fields are populated.
    """
    valid_by_chrom = _valid_by_chrom(valid_ranges)
    for solo in solos:
        # ---- primary path ----
        if solo.family and solo.family in family_to_sources:
            sources = family_to_sources[solo.family]
            resolved = _resolve_source_ervs_to_valid_ranges(sources, valid_ranges)
            if resolved:
                labels: set[str] = set()
                contributors: list[str] = []
                for r in resolved:
                    labels.update(r.probes)
                    if r.erv_id:
                        contributors.append(r.erv_id)
                if labels:
                    solo.probe_labels = sorted(labels)
                    solo.contributing_ervs = contributors
                    solo.label_source = "family"
                    continue
        # ---- fallback path ----
        nearest = _nearest_valid_erv(solo, valid_by_chrom, max_distance)
        if nearest and nearest.probes:
            solo.probe_labels = sorted(set(nearest.probes))
            solo.contributing_ervs = [nearest.erv_id] if nearest.erv_id else []
            solo.label_source = "nearest_erv"


# ------------------------------------------------------------------
# Writers
# ------------------------------------------------------------------
def write_solo_ltr_gff3(solos: list[SoloLTR], output_path: Path) -> None:
    """Write annotated solo LTRs as GFF3.

    Attributes per feature: ``ID``, ``family``, ``probe_labels`` (comma-
    separated), ``contributing_ervs`` (comma-separated), ``label_source``.
    An empty solos list still emits the header so Snakemake sees a
    valid file.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as fout:
        fout.write("##gff-version 3\n")
        fout.write("##source-version RetroSeek solo_ltr_integrator.py\n")
        for i, solo in enumerate(solos, start=1):
            attrs_parts = [
                f"ID=soloLTR_{i}",
                f"family={solo.family or 'unresolved'}",
                f"probe_labels={','.join(solo.probe_labels) if solo.probe_labels else 'none'}",
                f"contributing_ervs={','.join(solo.contributing_ervs) if solo.contributing_ervs else 'none'}",
                f"label_source={solo.label_source}",
            ]
            row = [
                solo.chrom,
                "LTR_retriever",
                "solo_LTR",
                str(solo.start),
                str(solo.end),
                ".",
                solo.strand,
                ".",
                ";".join(attrs_parts),
            ]
            fout.write("\t".join(row) + "\n")


def compute_solo_intact_ratio(
    solos: list[SoloLTR],
    valid_ranges: list[ValidRange],
    species: str,
) -> pd.DataFrame:
    """Return per-probe-family solo/intact counts + ratios.

    Two label modes are emitted:

    - ``exclusive`` — counts a solo/intact in exactly one family (its
      ``probes`` list has length 1).
    - ``shared``    — counts multi-labelled entries in every family they
      claim (length > 1).

    Rows with zero intact count get ``solo_to_intact_ratio = NaN``.
    """

    def _bump(
        bucket: dict[tuple[str, str], dict[str, int]],
        key: tuple[str, str],
        field_name: str,
        amount: int = 1,
    ) -> None:
        if key not in bucket:
            bucket[key] = {"solo_count": 0, "intact_count": 0}
        bucket[key][field_name] += amount

    bucket: dict[tuple[str, str], dict[str, int]] = {}
    for r in valid_ranges:
        if not r.probes:
            continue
        mode = "exclusive" if len(r.probes) == 1 else "shared"
        for p in r.probes:
            _bump(bucket, (p, mode), "intact_count")
    for solo in solos:
        if not solo.probe_labels:
            continue
        mode = "exclusive" if len(solo.probe_labels) == 1 else "shared"
        for p in solo.probe_labels:
            _bump(bucket, (p, mode), "solo_count")
    rows = []
    for (probe, mode), counts in sorted(bucket.items()):
        solo_count = counts["solo_count"]
        intact_count = counts["intact_count"]
        ratio = (solo_count / intact_count) if intact_count else float("nan")
        total = solo_count + intact_count
        solo_fraction = (solo_count / total) if total else 0.0
        rows.append(
            {
                "species": species,
                "probe_family": probe,
                "label_mode": mode,
                "intact_count": intact_count,
                "solo_count": solo_count,
                "total_integrations": total,
                "solo_to_intact_ratio": ratio,
                "solo_fraction": solo_fraction,
            }
        )
    return pd.DataFrame(rows)


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    """Entry point."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--nmtf-pass-list",
        type=Path,
        required=True,
        help="LTR_retriever .nmtf.pass.list (solo + truncated LTRs).",
    )
    parser.add_argument(
        "--ltr-library",
        type=Path,
        required=True,
        help="LTR_retriever .LTRlib.fa (consensus library with source-ERV headers).",
    )
    parser.add_argument(
        "--valid-ranges", type=Path, required=True, help="RetroSeek valid_ranges.gff3."
    )
    parser.add_argument(
        "--genome",
        required=True,
        help="Genome / species identifier (used in ratio CSV species column).",
    )
    parser.add_argument(
        "--output-gff3",
        type=Path,
        required=True,
        help="Output path for annotated solo-LTR GFF3.",
    )
    parser.add_argument(
        "--output-ratio-csv",
        type=Path,
        required=True,
        help="Output path for per-genome solo/intact ratio CSV.",
    )
    parser.add_argument(
        "--nearest-erv-max-distance",
        type=int,
        default=10000,
        help="bp window for the nearest-ERV fallback label propagation.",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    valid_ranges = parse_valid_ranges(args.valid_ranges)
    family_to_sources = parse_ltr_library_headers(args.ltr_library)
    solos = parse_nmtf_pass_list(args.nmtf_pass_list)

    logging.info(
        "Loaded %d valid ranges, %d consensus families, %d candidate solo LTRs",
        len(valid_ranges),
        len(family_to_sources),
        len(solos),
    )

    propagate_labels(
        solos=solos,
        valid_ranges=valid_ranges,
        family_to_sources=family_to_sources,
        max_distance=args.nearest_erv_max_distance,
    )

    n_family = sum(1 for s in solos if s.label_source == "family")
    n_nearest = sum(1 for s in solos if s.label_source == "nearest_erv")
    n_none = sum(1 for s in solos if s.label_source == "none")
    logging.info(
        "Label propagation — family: %d, nearest_erv: %d, unresolved: %d",
        n_family,
        n_nearest,
        n_none,
    )

    write_solo_ltr_gff3(solos, args.output_gff3)

    ratio_df = compute_solo_intact_ratio(solos, valid_ranges, species=args.genome)
    args.output_ratio_csv.parent.mkdir(parents=True, exist_ok=True)
    ratio_df.to_csv(args.output_ratio_csv, index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main())
