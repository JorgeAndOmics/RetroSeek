"""Probe extractor: fetch the probe sequences from NCBI and pickle them.

This script parses a CSV table containing probe gene metadata, retrieves their
corresponding sequences from an online database (e.g., NCBI GenBank), and serializes
the resulting data as pickled Python objects.

Modules:
    - `log`: the console line contract and job log (ADR-021).
    - `RetroSeek_class`: Contains the `RetroSeek` class representing a probe with metadata.
    - `seq_utils`: Provides functions for online sequence retrieval.
    - `utils`: Contains general-purpose utilities (e.g., string generator, pickling).
    - `defaults`: Stores configuration constants like directory paths and flags.

Main Workflow:
    1. Reads a probe metadata, comma-separated CSV file (columns: Label, Name, Abbreviation, Probe, Accession).
    2. Constructs a dictionary of `RetroSeek` instances keyed by accession ID.
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

import logging
from pathlib import Path

import pandas as pd

import defaults
import seq_utils
import utils
from log import OK, PipelineError, job_logging, run_main
from RetroSeeker_class import RetroSeeker

logger = logging.getLogger(__name__)


# =============================================================================
# 1. CSV Table Parser
# =============================================================================
def table_parser(input_csv_file: str | Path) -> dict[str, RetroSeeker]:
    """Parses a CSV file containing probe gene metadata and returns a dictionary of Object instances.

        Parameters
        ----------
            :param input_csv_file: Path to the input CSV file containing probe metadata.
                                   Expected columns: Label, Name, Abbreviation, Probe, Accession.

    Returns:
        -------
            :returns: probe_dict: A dictionary mapping accession IDs to `RetroSeeker` instances that encapsulate probe metadata.

    Raises:
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
            label=str(row["Label"]),
            virus=str(row["Name"]),
            abbreviation=str(row["Abbreviation"]),
            probe=str(row["Probe"]).upper(),
            accession=str(row["Accession"]),
            identifier=utils.random_string_generator(6),  # Random 6-char unique ID
        )
        for _, row in probe_table.iterrows()
    }

    return probe_dict


# =============================================================================
# 2. Main Execution: Extraction and Serialization
# =============================================================================
def main() -> None:
    """Fetch every probe named in the probe CSV and pickle them."""
    # Parse probe metadata from table
    probe_dict: dict[str, RetroSeeker] = table_parser(input_csv_file=defaults.PROBE_CSV)

    # Fetch sequence data from GenBank (or similar)
    probe_extraction: dict[str, RetroSeeker] | None = seq_utils.gb_executor(
        object_dict=probe_dict,
        online_database="protein",
    )

    if not probe_extraction:
        raise PipelineError(
            "no probe could be fetched from NCBI",
            hint="check the network, execution.entrez_email and the probe CSV accessions",
        )
    # A probe without its record would be skipped by every later search, so the
    # study would quietly lose it. Stop before anything is written instead.
    unfetched = [
        f"{probe.accession} ({probe.probe})"
        for probe in probe_extraction.values()
        if not probe.genbank
    ]
    if unfetched:
        raise PipelineError(
            f"{len(unfetched)} of {len(probe_extraction)} probes could not be fetched "
            f"from NCBI: {', '.join(unfetched[:10])}",
            hint="check these accessions in the probe CSV; if NCBI was only busy, "
            "run --probe-extractor again",
        )

    # Save extracted probes to serialized file
    utils.pickler(
        data=probe_extraction,
        output_directory_path=defaults.PATH_DICT["PICKLE_DIR"],
        output_file_name="probe_dict.pkl",
    )
    logger.log(OK, "%s probes fetched", f"{len(probe_extraction):,}")


if __name__ == "__main__":
    # The rule is a heavy one and must stay byte-identical, so it passes no log
    # path; the job log follows the LOG_DIR/<step>/all.log layout anyway.
    job_logging(
        Path(defaults.PATH_DICT["LOG_DIR"]) / "probe_extractor" / "all.log",
        "probe_extractor",
    )
    run_main(main)
