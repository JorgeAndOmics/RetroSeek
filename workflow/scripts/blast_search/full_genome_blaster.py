"""
Full Genome BLAST Executor
==========================

This script runs a tBLASTn search against a specified genome using a pre-serialized
dictionary of probe Objects. Results are serialized for downstream usage.

Modules:
    - `seq_utils`: Provides functions for executing sequence-related tasks.
    - `utils`: Contains helper functions.
    - `defaults`: Stores constant paths and parsed configuration values.
    - `log`: the console line contract and job log (ADR-021).

Main Workflow:
    1. Load pre-extracted probe dictionary from pickle.
    2. Run tBLASTn using those probes against a user-specified genome.
    3. Serialize the results as a `.pkl` file for reuse.

Requirements:
    - `probe_dict.pkl` must exist in `defaults.PICKLE_DIR`.
    - Genome database must exist in `defaults.SPECIES_DB`.

Usage:
    Run from the command line:
        python full_genome_blaster.py --genome MyGenomeName [--num_threads 2]
"""

# =============================================================================
# Imports and Logging Setup
# =============================================================================
import argparse
import logging
from pathlib import Path

import defaults
import seq_utils
import utils
from log import OK, job_logging, run_main
from RetroSeeker_class import RetroSeeker

logger = logging.getLogger(__name__)


# =============================================================================
# 1. Main Execution Block
# =============================================================================
def main(args: argparse.Namespace) -> None:
    """The work of one run; see the module docstring."""
    genome: str = args.genome
    num_threads: int = args.num_threads

    # -------------------------------------------------------------------------
    # 1.2 Load Probe Dictionary
    # -------------------------------------------------------------------------
    probe_dict: dict[str, RetroSeeker] = utils.unpickler(
        input_directory_path=defaults.PATH_DICT["PICKLE_DIR"],
        input_file_name="probe_dict.pkl",
    )

    # -------------------------------------------------------------------------
    # 1.3 Run tBLASTn Search
    # -------------------------------------------------------------------------
    tblastn_results: dict[str, RetroSeeker] | None = seq_utils.blast_retriever(
        object_dict=probe_dict,
        command="tblastn",
        genome=genome,
        input_database_path=defaults.PATH_DICT["SPECIES_DB"],
        num_threads=num_threads,
    )

    # -------------------------------------------------------------------------
    # 1.4 Serialize Results
    # -------------------------------------------------------------------------
    # A genome with zero tBLASTn hits yields None; persist an empty dict instead
    # so obj2dict.py loads every genome uniformly (calling .values() on None
    # would raise). The blast_pkl2parquet checkpoint then simply omits this
    # genome from species_with_hits() downstream.
    if tblastn_results is None:
        tblastn_results = {}

    utils.pickler(
        data=tblastn_results,
        output_directory_path=defaults.PATH_DICT["TBLASTN_PICKLE_DIR"],
        output_file_name=f"{genome}.pkl",
    )
    logger.log(OK, "%s tBLASTn hits", f"{len(tblastn_results):,}")


if __name__ == "__main__":
    # -------------------------------------------------------------------------
    # 1.1 Argument Parsing
    # -------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description="Performs tBLASTn on a genome.")
    parser.add_argument(
        "--genome",
        type=str,
        required=True,
        help="The name of the genome to perform tBLASTn on.",
    )
    parser.add_argument(
        "--num_threads",
        type=int,
        default=defaults.MAX_THREADPOOL_WORKERS,
        help="The number of threads to use for tBLASTn (default from config).",
    )
    args = parser.parse_args()
    # The rule is a heavy one and must stay byte-identical, so it passes no log
    # path; the job log follows the LOG_DIR/<step>/<genome>.log layout anyway.
    job_log = (
        Path(defaults.PATH_DICT["LOG_DIR"])
        / "full_genome_blaster"
        / f"{args.genome}.log"
    )
    job_logging(job_log, "full_genome_blaster")
    run_main(lambda: main(args))
