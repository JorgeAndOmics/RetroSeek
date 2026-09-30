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
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO

from Bio import SeqIO

from external import run_tool
from log import OK, PipelineError, job_logging, run_main
from tabular import write_csv

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


def _header_field(line: bytes) -> tuple[str, str] | None:
    """(field, value) for the header lines the labels use, else None.

    NAME is the HMMER tag; the RepeatMasker class is in Dfam's comment lines,
    `CC         Type: LTR` and `CC         SubType: ERVK`. Lines arrive as bytes
    and only these are decoded: most of the 11 GB are model rows.
    """
    if line.startswith(b"NAME "):
        return "NAME", line[5:].decode().strip()
    if line.startswith(b"CC ") and b":" in line:
        key, value = line[2:].decode(errors="replace").split(":", 1)
        if key.strip() in ("Type", "SubType"):
            return key.strip(), value.strip()
    return None


def _headers(handle: BinaryIO) -> Iterator[dict[str, str]]:
    """Each model's header fields in turn; a header ends at its `HMM` line."""
    fields: dict[str, str] = {}
    for line in handle:
        if line.startswith(b"HMM "):
            yield fields
            fields = {}
        elif field := _header_field(line):
            fields[field[0]] = field[1]


def read_model_classes(hmm: Path, wanted: set[str]) -> dict[str, str]:
    """The class of each `wanted` model: "LTR/ERVK", or "LINE" alone.

    Streams the file (about 11 GB) and stops once every wanted model is found.
    """
    classes: dict[str, str] = {}
    if not wanted:
        return classes
    with hmm.open("rb") as handle:
        for fields in _headers(handle):
            if fields.get("NAME") not in wanted:
                continue
            kinds = (fields.get("Type"), fields.get("SubType"))
            classes[fields["NAME"]] = "/".join(k for k in kinds if k)
            if len(classes) == len(wanted):
                break
    return classes


def _better(hit: Hit, current: Hit | None) -> bool:
    """Whether `hit` beats `current`: the higher score wins.

    A tie goes to the model name first in byte order, so the answer does not
    depend on the order the tables are read in.
    """
    if current is None:
        return True
    return (-hit.score, hit.model) < (-current.score, current.model)


def _table_hits(tblout: Path) -> Iterator[tuple[str, Hit]]:
    """(representative, hit) for every row of one `nhmmer --tblout` table."""
    for line in tblout.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        ali_from, ali_to, length = int(f[6]), int(f[7]), int(f[10])
        coverage = (abs(ali_to - ali_from) + 1) / length
        yield f[0], Hit(f[2], f[3], f[12], float(f[13]), coverage)


def best_hits(tblouts: list[Path]) -> dict[str, Hit]:
    """Each representative's highest-scoring model, over `nhmmer --tblout` tables."""
    best: dict[str, Hit] = {}
    for tblout in tblouts:
        for representative, hit in _table_hits(tblout):
            if _better(hit, best.get(representative)):
                best[representative] = hit
    return best


def label_rows(
    families: list[str], best: dict[str, Hit], classes: dict[str, str], release: str
) -> list[dict[str, object]]:
    """One row per family, blank where no curated model matched.

    Args:
        families: Qualified family names, `genome|family`.
        best: Each qualified name's best hit.
        classes: Each matched model's class.
        release: The Dfam release, recorded on every row.
    """
    rows: list[dict[str, object]] = []
    for family in families:
        genome, name = family.split("|", 1)
        hit = best.get(family)
        if hit is None:
            blank: dict[str, object] = dict.fromkeys(LABEL_COLUMNS, "")
            rows.append(
                blank | {"genome": genome, "ltr_family": name, "dfam_release": release}
            )
            continue
        rows.append(
            {
                "genome": genome,
                "ltr_family": name,
                "dfam_name": hit.model,
                "dfam_accession": hit.accession,
                "dfam_class": classes.get(hit.model, ""),
                "dfam_evalue": hit.evalue,
                "dfam_score": hit.score,
                "dfam_coverage": round(hit.coverage, 4),
                "dfam_release": release,
            }
        )
    return rows


def _read_fasta(fna: Path, wanted: set[str]) -> dict[str, str]:
    """The sequences of the `wanted` records."""
    # Biopython ships no type hints; validator.py silences the same call.
    records = SeqIO.parse(str(fna), "fasta")  # type: ignore[no-untyped-call]
    return {record.id: str(record.seq) for record in records if record.id in wanted}


def _representatives(summary_csv: Path) -> dict[str, str]:
    """Each family's representative arm -> the family's name (one row per family)."""
    with summary_csv.open(newline="", encoding="utf-8") as table:
        return {
            row["representative"]: row["ltr_family"] for row in csv.DictReader(table)
        }


def write_representatives(pairs: list[tuple[str, Path, Path]], out: Path) -> list[str]:
    """Write each family's representative arm under `genome|family`.

    Args:
        pairs: For each genome, its name, its family summary and its bait FASTA.
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
        for genome, summary_csv, bait_fna in pairs:
            reps = _representatives(summary_csv)
            sequences = _read_fasta(bait_fna, set(reps))
            missing = sorted(reps[arm] for arm in set(reps) - set(sequences))
            if missing:
                raise PipelineError(
                    f"{genome}: the representatives of {', '.join(missing)} are not "
                    f"in {bait_fna.name}",
                    hint="the family summary and the bait must come from one run",
                )
            for arm, sequence in sequences.items():
                handle.write(f">{genome}|{reps[arm]}\n{sequence}\n")
                names.append(f"{genome}|{reps[arm]}")
    return names


def split_models(hmm: Path, parts: int, workdir: Path) -> list[Path]:
    """Deal the models of `hmm` in turn into `parts` files; a model ends at `//`.

    Returns only the files that got a model: nhmmer refuses an empty one.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    paths = [workdir / f"models.{i}.hmm" for i in range(parts)]
    handles = [path.open("wb") for path in paths]
    try:
        turn = 0
        with hmm.open("rb") as source:
            for line in source:
                handles[turn].write(line)
                if line.startswith(b"//"):
                    turn = (turn + 1) % parts
    finally:
        for handle in handles:
            handle.close()
    for path in paths:
        if path.stat().st_size == 0:
            path.unlink()
    return [path for path in paths if path.exists()]


def _nhmmer(models: Path, representatives: Path, tblout: Path) -> None:
    """One single-threaded nhmmer: every model in `models` against the representatives."""
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
            "1",
            str(models),
            str(representatives),
        ]
    )


def run_nhmmer(
    hmm: Path, representatives: Path, workdir: Path, threads: int
) -> list[Path]:
    """Search the representatives with every model; returns the `--tblout` tables.

    nhmmer's own threads share out the target sequences, and a few hundred short
    representatives make a single block: on the model genomes 8 threads kept 1.5
    cores busy, for 2.6 hours. So the models are dealt into `threads` chunks in
    `workdir` (a copy of the 11 GB file, removed afterwards) and each chunk gets
    its own single-threaded nhmmer.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    chunks = split_models(hmm, threads, workdir) if threads > 1 else [hmm]
    tblouts = [workdir / f"dfam.{i}.tbl" for i in range(len(chunks))]
    try:
        with ThreadPoolExecutor(max_workers=len(chunks)) as pool:
            list(pool.map(_nhmmer, chunks, [representatives] * len(chunks), tblouts))
    finally:
        for chunk in chunks:
            if chunk != hmm:
                chunk.unlink(missing_ok=True)
    return tblouts


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
        (
            g,
            args.table_dir / f"{g}.ltr_family_summary.csv",
            args.bait_dir / f"{g}.bait.fna",
        )
        for g in sorted(args.genomes)
    ]
    reps = args.workdir / "representatives.fna"
    families = write_representatives(pairs, reps)
    best = best_hits(run_nhmmer(args.dfam_hmm, reps, args.workdir, args.threads))
    classes = read_model_classes(args.dfam_hmm, {hit.model for hit in best.values()})
    rows = label_rows(families, best, classes, args.release)
    write_csv(rows, LABEL_COLUMNS, args.out)
    logger.log(
        OK,
        "%s of %s LTR families matched a curated Dfam model",
        f"{len(best):,}",
        f"{len(families):,}",
    )


if __name__ == "__main__":
    run_main(main)
