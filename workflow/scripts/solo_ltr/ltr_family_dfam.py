# =============================================================================
# ltr_family_dfam.py: label LTR families with their best Dfam model (ADR-023)
# =============================================================================
# Optional, off by default (`solo_ltr.families.dfam`). Each family's representative
# arm is searched with every curated Dfam model (`nhmmer`, the models as queries,
# Dfam's own gathering thresholds), and the best-scoring model gives the family an
# outside name and class, e.g. IAPLTR1_Mm, LTR/ERVK. The label is evidence to
# read, never a filter: RetroSeek still finds every element and solo itself.
#
# Coverage is counted over the representative, since `--tblout` has no model
# length. Dfam's curated models cover human and mouse well and the bats mostly
# through families shared across mammals, so a blank label means "no curated
# model matched", not "not an ERV".
# =============================================================================

"""Label every genome's LTR families with the best-matching curated Dfam model."""

from __future__ import annotations

import argparse
import csv
import logging
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from ltr_families import write_csv

from external import run_tool
from log import OK, job_logging, run_main

logger = logging.getLogger(__name__)

NHMMER = "nhmmer"
LABEL_COLUMNS = [
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


def _ct_value(field: str) -> str:
    """The value of a Dfam `CT` line: `Type; LTR;` gives `LTR`."""
    return field.split(";", 1)[-1].strip().rstrip(";").strip()


def read_model_headers(hmm: Path) -> dict[str, tuple[str, str]]:
    """Each model's (accession, class) from a Dfam HMM file, read as a stream.

    The class joins the model's `CT` values: `Type; LTR` and `SubType; ERVK` give
    `LTR/ERVK`. HMMER header fields are tags left-justified in six columns.
    """
    models: dict[str, tuple[str, str]] = {}
    name, accession = "", ""
    classes: list[str] = []
    with hmm.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            tag = line[:6].strip()
            if tag == "NAME":
                name, accession, classes = line[6:].strip(), "", []
            elif tag == "ACC":
                accession = line[6:].strip()
            elif tag == "CT":
                classes.append(_ct_value(line[6:]))
            elif tag == "HMM" and name:
                models[name] = (accession, "/".join(classes))
                name = ""
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
    """One row per family, blank where no curated model matched."""
    rows: list[dict[str, object]] = []
    for family in families:
        hit = best.get(family)
        if hit is None:
            blank: dict[str, object] = dict.fromkeys(LABEL_COLUMNS, "")
            rows.append(blank | {"ltr_family": family})
            continue
        rows.append(
            {
                "ltr_family": family,
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


def write_representatives(pairs: list[tuple[Path, Path]], out: Path) -> list[str]:
    """Write each family's representative arm under the family's name.

    Args:
        pairs: For each genome, its family table and its bait FASTA.
        out: The FASTA to write.

    Returns:
        The family names written, in order.
    """
    names = []
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        for families_csv, bait_fna in pairs:
            with families_csv.open(newline="", encoding="utf-8") as table:
                reps = {
                    row["arm"]: row["ltr_family"]
                    for row in csv.DictReader(table)
                    if row["representative"] == "True"
                }
            for arm, sequence in _read_fasta(bait_fna, set(reps)).items():
                handle.write(f">{reps[arm]}\n{sequence}\n")
                names.append(reps[arm])
    return names


def run_nhmmer(hmm: Path, representatives: Path, workdir: Path, threads: int) -> Path:
    """Search the representatives with every model at its gathering threshold."""
    workdir.mkdir(parents=True, exist_ok=True)
    tblout = workdir / "dfam.tbl"
    run_tool(
        [
            NHMMER,
            "--cut_ga",
            "--noali",
            "--tblout",
            str(tblout),
            "-o",
            str(workdir / "nhmmer.out"),
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
        (args.table_dir / f"{g}.ltr_families.csv", args.bait_dir / f"{g}.bait.fna")
        for g in sorted(args.genomes)
    ]
    reps = args.workdir / "representatives.fna"
    families = write_representatives(pairs, reps)
    best = best_hits(run_nhmmer(args.dfam_hmm, reps, args.workdir, args.threads))
    rows = label_rows(families, best, read_model_headers(args.dfam_hmm))
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
