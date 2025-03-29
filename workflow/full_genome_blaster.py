"""
Full Genome BLAST Executor
==========================

This script runs a tBLASTn search against a specified genome using a pre-serialized
dictionary of probe Objects. Results are serialized for downstream usage.

Modules:
    - `seq_utils`: Provides functions for executing sequence-related tasks.
    - `utils`: Contains helper functions.
    - `defaults`: Stores constant paths and parsed configuration values.
    - `colored_logging`: Enables structured, configurable, color-coded logging output.

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

from colored_logging import colored_logging
import seq_utils
import utils
import defaults

# =============================================================================
# 1. Main Execution Block
# =============================================================================
if __name__ == '__main__':
    # Initialize logging
    colored_logging(log_file_name='full_genome_blaster.txt')

    # -------------------------------------------------------------------------
    # 1.1 Argument Parsing
    # -------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Performs tBLASTn on a genome.')
    parser.add_argument(
        '--genome',
        type=str,
        required=True,
        help='The name of the genome to perform tBLASTn on.'
    )
    parser.add_argument(
        '--num_threads',
        type=int,
        default=defaults.MAX_THREADPOOL_WORKERS,
        help='The number of threads to use for tBLASTn (default from config).'
    )
    args = parser.parse_args()

    genome: str = args.genome
    num_threads: int = args.num_threads

    # -------------------------------------------------------------------------
    # 1.2 Load Probe Dictionary
    # -------------------------------------------------------------------------
    probe_dict: dict[str, object] = utils.unpickler(
        input_directory_path=defaults.PICKLE_DIR,
        input_file_name='probe_dict.pkl'
    )

    # -------------------------------------------------------------------------
    # 1.3 Run tBLASTn Search
    # -------------------------------------------------------------------------
    tblastn_results: dict[str, object] = seq_utils.blast_retriever(
        object_dict=probe_dict,
        command='tblastn',
        genome=genome,
        input_database_path=defaults.SPECIES_DB,
        num_threads=num_threads,
        display_full_info=False
    )

    # -------------------------------------------------------------------------
    # 1.4 Serialize Results
    # -------------------------------------------------------------------------
    utils.pickler(
        data=tblastn_results,
        output_directory_path=defaults.TBLASTN_PICKLE_DIR,
        output_file_name=f'{genome}.pkl'
    )
