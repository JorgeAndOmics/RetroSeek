"""
Probe Extractor Script
======================

This script parses a CSV table containing probe gene metadata, retrieves their
corresponding sequences from an online database (e.g., NCBI GenBank), and serializes
the resulting data as pickled Python objects.

Modules:
    - `colored_logging`: Custom logging with color-coded console/file output.
    - `object_class`: Contains the `Object` class representing a probe with metadata.
    - `seq_utils`: Provides functions for online sequence retrieval.
    - `utils`: Contains general-purpose utilities (e.g., string generator, pickling).
    - `defaults`: Stores configuration constants like directory paths and flags.

Main Workflow:
    1. Reads a probe metadata, comma-separated CSV file (columns: Family, Name, Abbreviation, Probe, Accession).
    2. Constructs a dictionary of `Object` instances keyed by accession ID.
    3. Fetches sequences from a remote database using the given metadata.
    4. Serializes the extracted objects to a pickle file for downstream use.

Requirements:
    - The input CSV must be located in `defaults.TABLE_INPUT_DIR`.
    - Output pickled data will be saved to `defaults.PICKLE_DIR`.

Usage:
    Run as a standalone script:
        python probe_extractor.py

    Ensure the required CSV (probes.csv) exists at the designated location (config.yaml).
"""

# =============================================================================
# Imports and Logging Setup
# =============================================================================
from colored_logging import colored_logging
import pandas as pd
import logging
import os

from RetroSeeker_class import RetroSeeker
import seq_utils
import utils
import defaults

# =============================================================================
# 1. CSV Table Parser
# =============================================================================
def table_parser(input_csv_file: str) -> dict[str, RetroSeeker]:
    """
    Parses a CSV file containing probe gene metadata and returns a dictionary of Object instances.

        Parameters
        ----------
            :param input_csv_file: Path to the input CSV file containing probe metadata.
                                   Expected columns: Family, Name, Abbreviation, Probe, Accession.

        Returns
        -------
            :returns: probe_dict: A dictionary mapping accession IDs to `Object` instances that encapsulate probe metadata.

        Raises
        ------
            :raises FileNotFoundError: If the specified CSV file does not exist.
            :raises pd.errors.ParserError: If the CSV content is malformed.
            :raises KeyError: If required columns are missing from the CSV.
    """
    # Read the CSV file into a DataFrame
    probe_table = pd.read_csv(input_csv_file)

    # Build dictionary mapping accession IDs to Object instances
    probe_dict: dict[str, RetroSeeker] = {
        str(row["Accession"]): RetroSeeker(
            family=str(row['Family']),
            virus=str(row['Name']),
            abbreviation=str(row['Abbreviation']),
            probe=str(row['Probe']),
            accession=str(row['Accession']),
            identifier=utils.random_string_generator(6)  # Random 6-char unique ID
        )
        for _, row in probe_table.iterrows()
    }

    return probe_dict


# =============================================================================
# 2. Main Execution: Extraction and Serialization
# =============================================================================
if __name__ == '__main__':
    # Initialize colored logging and write to log
    colored_logging(log_file_name='probe_extractor.txt')

    # Parse probe metadata from table
    probe_dict: dict[str, RetroSeeker] = table_parser(
        input_csv_file=defaults.PROBE_CSV)

    # Fetch sequence data from GenBank (or similar)
    probe_extraction: dict[str, RetroSeeker] = seq_utils.gb_executor(
        object_dict=probe_dict,
        online_database='protein',
        display_full_info=defaults.DISPLAY_OPERATION_INFO,
        display_warning=defaults.DISPLAY_REQUESTS_WARNING
    )

    # Save extracted probes to serialized file
    utils.pickler(
        data=probe_extraction,
        output_directory_path=defaults.PATH_DICT['PICKLE_DIR'],
        output_file_name='probe_dict.pkl'
    )
