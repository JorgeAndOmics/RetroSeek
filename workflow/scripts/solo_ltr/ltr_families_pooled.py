# =============================================================================
# ltr_families_pooled.py: LTR families across every genome at once (ADR-023)
# =============================================================================
# The same clustering as ltr_families.py, over all genomes' bait arms together,
# so a family shared by two bats is one family (`Pool_F001`...). Arm names are
# prefixed with their genome, because chromosome names can repeat between
# assemblies. Genomes are pooled in sorted order: cd-hit-est's result can depend
# on input order, and the config's order must not change the families.
# =============================================================================

"""Group every genome's solo-LTR bait arms into pooled LTR families."""

from __future__ import annotations

import argparse
import logging
from collections import Counter
from pathlib import Path
from typing import TextIO

from bait_builder import Arm
from ltr_families import Family, name_families, read_bait_bed, run_cdhit

from log import OK, job_logging, run_main
from tabular import write_csv

logger = logging.getLogger(__name__)

POOL_PREFIX = "Pool"
POOLED_COLUMNS = ["genome", "arm", "pool_family", "representative"]
POOLED_SUMMARY_COLUMNS = [
    "pool_family",
    "n_arms",
    "n_genomes",
    "arms_per_genome",
    "representative",
]


def _copy_prefixed(fna: Path, genome: str, out: TextIO) -> None:
    """Copy one genome's bait FASTA into `out`, each name prefixed `genome|`."""
    with fna.open(encoding="utf-8") as bait:
        for line in bait:
            out.write(f">{genome}|{line[1:]}" if line.startswith(">") else line)


def pool_bait(bait_dir: Path, genomes: list[str], out_fna: Path) -> dict[str, Arm]:
    """Write every genome's bait arms into one FASTA, names prefixed by genome.

    Returns the pooled arms by prefixed name, placed on `genome|seqname` so that
    family order stays deterministic.
    """
    arms: dict[str, Arm] = {}
    out_fna.parent.mkdir(parents=True, exist_ok=True)
    with out_fna.open("w", encoding="utf-8") as out:
        for genome in sorted(genomes):
            _copy_prefixed(bait_dir / f"{genome}.bait.fna", genome, out)
            for name, arm in read_bait_bed(bait_dir / f"{genome}.bait.bed").items():
                pooled = f"{genome}|{name}"
                seqname = f"{genome}|{arm.seqname}"
                arms[pooled] = Arm(seqname, arm.start, arm.end, arm.parent, arm.arm)
    return arms


def pooled_rows(families: list[Family]) -> list[dict[str, object]]:
    """One row per arm: its genome, its own name, and its pooled family."""
    rows: list[dict[str, object]] = []
    for family in families:
        for member in family.members:
            genome, arm = member.split("|", 1)
            rows.append(
                {
                    "genome": genome,
                    "arm": arm,
                    "pool_family": family.name,
                    "representative": member == family.representative,
                }
            )
    return rows


def pooled_summary(families: list[Family]) -> list[dict[str, object]]:
    """One row per pooled family: its size and which genomes it spans."""
    rows: list[dict[str, object]] = []
    for family in families:
        per_genome = Counter(member.split("|", 1)[0] for member in family.members)
        rows.append(
            {
                "pool_family": family.name,
                "n_arms": len(family.members),
                "n_genomes": len(per_genome),
                "arms_per_genome": ";".join(
                    f"{g}:{n}" for g, n in sorted(per_genome.items())
                ),
                "representative": family.representative,
            }
        )
    return rows


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    """The command line: where the bait is, which genomes, and the outputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bait-dir", type=Path, required=True)
    parser.add_argument("--genomes", nargs="+", required=True)
    parser.add_argument("--identity", type=float, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--out-families", type=Path, required=True)
    parser.add_argument("--out-summary", type=Path, required=True)
    parser.add_argument("--log", type=Path, help="job log (the Snakemake log: path)")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    """Command-line entry: pool every genome's bait arms and write the families."""
    args = _parse_args(argv)
    job_logging(args.log, "solo_family_pooled")
    pooled_fna = args.workdir / "pooled.bait.fna"
    arms = pool_bait(args.bait_dir, args.genomes, pooled_fna)
    clusters = run_cdhit(pooled_fna, args.workdir, args.identity, args.threads)
    families = name_families(clusters, arms, POOL_PREFIX)
    write_csv(pooled_rows(families), POOLED_COLUMNS, args.out_families)
    summary = pooled_summary(families)
    write_csv(summary, POOLED_SUMMARY_COLUMNS, args.out_summary)
    shared = sum(int(row["n_genomes"] != 1) for row in summary)
    logger.log(
        OK,
        "%s arms from %s genomes in %s pooled families, %s shared by two or more",
        f"{len(arms):,}",
        len(args.genomes),
        f"{len(families):,}",
        f"{shared:,}",
    )


if __name__ == "__main__":
    run_main(main)
