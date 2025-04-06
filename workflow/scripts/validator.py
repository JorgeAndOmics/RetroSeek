# -------------------
# DEPENDENCIES
# -------------------

import os
import re
import sys
import time
import argparse
import subprocess
import logging
from typing import List, Optional

import yamale
import pandas as pd
import pandera as pa
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import defaults
import validator
from colored_logging import colored_logging


# -----------------------------
# YAML VALIDATION
# -----------------------------

def yaml_validator(yaml_file: str, yaml_schema: str) -> bool:
    """
    Validates the YAML configuration file against a Yamale schema.

    Parameters
    ----------
    yaml_file : str
        Path to the YAML configuration file.
    yaml_schema : str
        Path to the YAML schema file.

    Returns
    -------
    bool
        True if the YAML is valid, False otherwise.
    """
    try:
        schema = yamale.make_schema(yaml_schema)
        data = yamale.make_data(yaml_file)
        yamale.validate(schema, data)
        logging.info('YAML configuration file is valid.')
        return True

    except yamale.YamaleError as e:
        for results in e.results:
            for error in results.errors:
                logging.warning(f"Configuration error: {error}")
        return False


# -----------------------------
# CSV VALIDATION
# -----------------------------

def csv_validator(csv_file: str) -> bool:
    """
    Validates the CSV file for required fields and checks NCBI accession validity.

    Parameters
    ----------
    csv_file : str
        Path to the CSV file.

    Returns
    -------
    bool
        True if CSV passes validation, False otherwise.
    """

    def check_ncbi(series: pd.Series) -> pd.Series:
        """Check if accession IDs exist in the NCBI protein database."""
        logging.debug('Checking NCBI entries...')
        def ncbi_exists(acc: str) -> bool:
            time.sleep(0.3)
            Entrez.email = defaults.ENTREZ_EMAIL
            try:
                with Entrez.esearch(db='protein', term=acc) as handle:
                    record = Entrez.read(handle)
                    return int(record["Count"]) > 0
            except Exception:
                return False

        return series.apply(ncbi_exists)

    try:
        df = pd.read_csv(csv_file)

        schema = pa.DataFrameSchema({
            'Label': pa.Column(str, nullable=False),
            'Name': pa.Column(str, nullable=False),
            'Abbreviation': pa.Column(str, nullable=False),
            'Probe': pa.Column(str, nullable=False),
            'Accession': pa.Column(
                str, nullable=False,
                checks=pa.Check(check_ncbi, element_wise=False)
            )
        })

        schema.validate(df)
        logging.info('CSV input file is valid.')
        return True

    except (pa.errors.SchemaError, FileNotFoundError, KeyError) as e:
        logging.warning(f'CSV input error: {e}')
        return False


# -----------------------------
# FASTA VALIDATION
# -----------------------------

def fasta_validator(fasta_file: str) -> bool:
    """
    Validates a FASTA file for content and headers.

    Parameters
    ----------
    fasta_file : str
        Path to the FASTA file.

    Returns
    -------
    bool
        True if valid, False otherwise.
    """
    if not os.path.exists(fasta_file):
        logging.warning(f"FASTA file does not exist: {fasta_file}")
        return False

    try:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        if not records:
            logging.warning("FASTA file is empty or has no valid records.")
            return False

        for i, record in enumerate(records):
            if not record.id:
                logging.warning(f"Record {i + 1} is missing a header.")
                return False

        logging.info(f'FASTA file {fasta_file} is valid.')
        return True

    except Exception as e:
        logging.warning(f"Error parsing FASTA: {e}")
        return False


# -----------------------------
# TOOL AVAILABILITY CHECK
# -----------------------------

def validate_programs() -> bool:
    """
    Validates required external programs are available in the system path.

    Returns
    -------
    bool
        True if all required programs are available, False otherwise.
    """

    def check_version(cmd: str, version_cmd: str) -> bool:
        try:
            subprocess.run([cmd, version_cmd], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False

    checks = {
        'BLAST+': check_version('tblastn', '-version'),
        'Genometools': check_version('gt', '--version'),
        'Datasets': check_version('datasets', '--version')
    }

    for tool, status in checks.items():
        log_fn = logging.info if status else logging.warning
        log_fn(f"{tool} {'is installed.' if status else 'not found. Please install it.'}")

    return all(checks.values())


# -----------------------------
# NCBI API KEY VALIDATION
# -----------------------------

def validate_ncbi_key() -> None:
    """
    Ensures that an NCBI API key is available in the environment.
    Prompts the user to enter one if not present.
    """
    if 'NCBI_API_KEY' not in os.environ:
        logging.warning('No NCBI API key found. You can set it with: export NCBI_API_KEY="your_api_key".')
        if (
            api_key := input('Enter your NCBI API key [Leave empty to skip]: ')
            or None
        ):
            os.environ['NCBI_API_KEY'] = api_key
            logging.info('NCBI API key set in environment variables.')
        else:
            logging.warning('No NCBI API key provided. Expect slower NCBI retrievals.')
    else:
        logging.info("NCBI API key is set in the environment variables.")


# -----------------------------
# MASTER VALIDATION ENTRYPOINT
# -----------------------------

def main_validator(fasta_files: Optional[List[str]]) -> bool:
    """
    Orchestrates validation for YAML, CSV, FASTA, external tools, and API key.

    Parameters
    ----------
    fasta_files : list of str or None
        List of FASTA file paths to validate.

    Returns
    -------
    bool
        True if all checks pass, False otherwise.
    """
    logging.debug('Starting input validation process...')

    yaml_ok = yaml_validator(
        yaml_schema=os.path.join(defaults.PATH_DICT['CONFIG_DIR'], 'schema.yaml'),
        yaml_file=defaults.CONFIG_FILE
    )

    csv_ok = csv_validator(csv_file=defaults.config['input']['probe_csv'])

    if not defaults.USE_SPECIES_DICT and fasta_files:
        fasta_results = [fasta_validator(f) for f in fasta_files]
        fasta_ok = all(fasta_results)
    else:
        fasta_ok = True

    programs_ok = validate_programs()

    validate_ncbi_key()

    return all([yaml_ok, csv_ok, fasta_ok, programs_ok])


# -----------------------------
# INTERACTIVE CONFIRMATION
# -----------------------------

def green_light(all_valid: bool) -> bool:
    """
    Asks user to confirm whether to proceed if all validations passed.

    Parameters
    ----------
    all_valid : bool
        Whether all validations were successful.

    Returns
    -------
    bool
        True if user wants to proceed, False otherwise.
    """
    if not all_valid:
        logging.warning(f'Settings validation failed. Check logs at {defaults.PATH_DICT["LOG_DIR"]}.')
        return False

    logging.info('All systems green, ready to rock.')
    time.sleep(0.1)

    proceed = input('Proceed [Y/n]: ') or 'Y'
    if proceed.upper() == 'Y':
        logging.info('RetroSeek started. Depending on your system, this may take some time.')
        return True
    elif proceed.upper() == 'N':
        logging.warning('Workflow aborted by user.')
        return False
    else:
        logging.warning('Invalid input. Please try again.')
        green_light(all_valid)


# -----------------------------
# ENTRYPOINT
# -----------------------------

def validation_run(fasta_files: Optional[List[str]] = None) -> bool:
    """
    CLI entrypoint to trigger validation routines and prompt user to continue.

    Parameters
    ----------
    fasta_files : list of str, optional
        List of FASTA file paths to validate.

    Returns
    -------
    bool
        True if user confirms execution after passing validation, False otherwise.
    """
    colored_logging(log_file_name='validator.log')

    fasta_files_list = fasta_files or []
    all_valid = main_validator(fasta_files=fasta_files_list)

    return green_light(all_valid=all_valid)
