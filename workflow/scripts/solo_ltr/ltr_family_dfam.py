# =============================================================================
# ltr_family_dfam.py: label LTR families with their best Dfam model (ADR-023)
# =============================================================================
# Optional, off by default (`solo_ltr.families.dfam`). Each family's representative
# arm is searched with every curated Dfam model (`nhmmer`, the models as queries),
# and the best-scoring model gives the family an outside name and class, e.g.
# IAPLTR1_Mm, LTR/ERVK. The label is evidence to read, never a filter: RetroSeek
# still finds every element and solo itself.
#
# One E-value for every model, not Dfam's gathering thresholds: in Dfam 4.0 about
# 5,600 of the 30,646 curated models have no model-level GA line (`--cut_ga` then
# stops nhmmer), and where there is one it equals the strict TC, while the
# per-taxon thresholds sit in TH lines. At 1e-5 per model the whole search expects
# well under one chance hit.
#
# Coverage is counted over the representative, since `--tblout` has no model
# length. Families are named `genome|family` inside the search, because two
# genomes can share a family code (Canis_lupus_familiaris and Canis_lupus_dingo
# are both Clup).
# =============================================================================

"""Label every genome's LTR families with the best-matching curated Dfam model."""

from __future__ import annotations

import argparse
import csv
import logging
import os
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from Bio import SeqIO
from ltr_families import write_csv

from external import run_tool
from log import OK, PipelineError, job_logging, run_main

logger = logging.getLogger(__name__)

NHMMER = "nhmmer"
EVALUE = 1e-5  # per model; see the header for why not --cut_ga
LABEL_COLUMNS = [
    "genome",
    "ltr_family",
    "dfam_name",
    "dfam_accession",
    "dfam_class",
    "dfam_evalue",
    "dfam_score",
    "dfam_coverage",
    "dfam_release",
]


@dataclass(frozen=True)
class Hit:
    """The best Dfam model found in one representative arm."""

    model: str
    accession: str
    evalue: str
    score: float
    coverage: float  # share of the representative covered by the alignment


def _header_field(line: str) -> tuple[str, str] | None:
    """(field, value) for the header lines the labels use, else None.

    NAME and ACC are HMMER tags; the RepeatMasker class is in Dfam's comment
    lines, `CC         Type: LTR` and `CC         SubType: ERVK`.
    """
    if line.startswith(("NAME ", "ACC ")):
        tag, value = line.split(maxsplit=1)
        return tag, value.strip()
    if line.startswith("CC ") and ":" in line:
        key, value = line[2:].split(":", 1)
        if key.strip() in ("Type", "SubType"):
            return key.strip(), value.strip()
    return None


def _headers(handle: TextIO) -> Iterator[dict[str, str]]:
    """Each model's header fields in turn; a header ends at its `HMM` line."""
    fields: dict[str, str] = {}
    for line in handle:
        if line.startswith("HMM "):
            yield fields
            fields = {}
        elif field := _header_field(line):
            fields[field[0]] = field[1]


def _model_class(fields: dict[str, str]) -> str:
    """RepeatMasker class and subclass joined: "LTR/ERVK", or "LINE" alone."""
    return "/".join(c for c in (fields.get("Type"), fields.get("SubType")) if c)


def read_model_headers(hmm: Path, wanted: set[str]) -> dict[str, tuple[str, str]]:
    """The (accession, class) of each `wanted` model, e.g. ("DF000004162.1", "LTR/ERVK").

    Streams the file (about 11 GB) and stops once every wanted model is found.
    """
    models: dict[str, tuple[str, str]] = {}
    if not wanted:
        return models
    with hmm.open(encoding="utf-8", errors="replace") as handle:
        for fields in _headers(handle):
            if fields.get("NAME") not in wanted:
                continue
            models[fields["NAME"]] = (fields.get("ACC", ""), _model_class(fields))
            if len(models) == len(wanted):
                break
    return models


def best_hits(tblout: Path) -> dict[str, Hit]:
    """Each representative's highest-scoring model, from `nhmmer --tblout`."""
    best: dict[str, Hit] = {}
    for line in tblout.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        ali_from, ali_to, length = int(f[6]), int(f[7]), int(f[10])
        hit = Hit(
            f[2], f[3], f[12], float(f[13]), (abs(ali_to - ali_from) + 1) / length
        )
        if f[0] not in best or hit.score > best[f[0]].score:
            best[f[0]] = hit
    return best


def label_rows(
    families: list[str], best: dict[str, Hit], models: dict[str, tuple[str, str]]
) -> list[dict[str, object]]:
    """One row per family, blank where no curated model matched.

    Args:
        families: Qualified family names, `genome|family`.
        best: Each qualified name's best hit.
        models: Each matched model's (accession, class).
    """
    rows: list[dict[str, object]] = []
    for family in families:
        genome, name = family.split("|", 1)
        hit = best.get(family)
        if hit is None:
            blank: dict[str, object] = dict.fromkeys(LABEL_COLUMNS, "")
            rows.append(blank | {"genome": genome, "ltr_family": name})
            continue
        rows.append(
            {
                "genome": genome,
                "ltr_family": name,
                "dfam_name": hit.model,
                "dfam_accession": hit.accession,
                "dfam_class": models.get(hit.model, ("", ""))[1],
                "dfam_evalue": hit.evalue,
                "dfam_score": hit.score,
                "dfam_coverage": round(hit.coverage, 4),
            }
        )
    return rows


def _read_fasta(fna: Path, wanted: set[str]) -> dict[str, str]:
    """The sequences of the `wanted` records."""
    # Biopython ships no type hints; validator.py silences the same call.
    records = SeqIO.parse(str(fna), "fasta")  # type: ignore[no-untyped-call]
    return {record.id: str(record.seq) for record in records if record.id in wanted}


def _representatives(families_csv: Path) -> dict[str, str]:
    """Each family's representative arm -> the family's name."""
    with families_csv.open(newline="", encoding="utf-8") as table:
        return {
            row["arm"]: row["ltr_family"]
            for row in csv.DictReader(table)
            if row["representative"] == "True"
        }


def write_representatives(pairs: list[tuple[str, Path, Path]], out: Path) -> list[str]:
    """Write each family's representative arm under `genome|family`.

    Args:
        pairs: For each genome, its name, its family table and its bait FASTA.
        out: The FASTA to write.

    Returns:
        The qualified family names written, in order.

    Raises:
        PipelineError: If a representative is missing from its bait FASTA; the
            two files would then come from different runs.
    """
    names = []
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        for genome, families_csv, bait_fna in pairs:
            reps = _representatives(families_csv)
            sequences = _read_fasta(bait_fna, set(reps))
            missing = sorted(reps[arm] for arm in set(reps) - set(sequences))
            if missing:
                raise PipelineError(
                    f"{genome}: the representatives of {', '.join(missing)} are not "
                    f"in {bait_fna.name}",
                    hint="the family table and the bait must come from one run",
                )
            for arm, sequence in sequences.items():
                handle.write(f">{genome}|{reps[arm]}\n{sequence}\n")
                names.append(f"{genome}|{reps[arm]}")
    return names


def run_nhmmer(hmm: Path, representatives: Path, workdir: Path, threads: int) -> Path:
    """Search the representatives with every model; returns the `--tblout` table."""
    workdir.mkdir(parents=True, exist_ok=True)
    tblout = workdir / "dfam.tbl"
    run_tool(
        [
            NHMMER,
            "-E",
            str(EVALUE),
            "--noali",
            "--tblout",
            str(tblout),
            "-o",
            os.devnull,  # the per-model report is large and nothing reads it
            "--cpu",
            str(threads),
            str(hmm),
            str(representatives),
        ]
    )
    return tblout


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    """The command line: the genomes, where their tables are, Dfam and the output."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genomes", nargs="+", required=True)
    parser.add_argument("--table-dir", type=Path, required=True)
    parser.add_argument("--bait-dir", type=Path, required=True)
    parser.add_argument("--dfam-hmm", type=Path, required=True)
    parser.add_argument("--release", required=True, help="Dfam release, e.g. 4.0")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--log", type=Path, help="job log (the Snakemake log: path)")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    """Command-line entry: label every genome's families and write one table."""
    args = _parse_args(argv)
    job_logging(args.log, "solo_family_dfam")
    pairs = [
        (g, args.table_dir / f"{g}.ltr_families.csv", args.bait_dir / f"{g}.bait.fna")
        for g in sorted(args.genomes)
    ]
    reps = args.workdir / "representatives.fna"
    families = write_representatives(pairs, reps)
    best = best_hits(run_nhmmer(args.dfam_hmm, reps, args.workdir, args.threads))
    models = read_model_headers(args.dfam_hmm, {hit.model for hit in best.values()})
    rows = label_rows(families, best, models)
    for row in rows:
        row["dfam_release"] = args.release
    write_csv(rows, LABEL_COLUMNS, args.out)
    logger.log(
        OK,
        "%s of %s LTR families matched a curated Dfam model",
        f"{len(best):,}",
        f"{len(families):,}",
    )


if __name__ == "__main__":
    run_main(main)
