import defaults

import pandera as pa
import pandas as pd
import subprocess
import argparse
import yamale
import time
import os

from Bio.SeqRecord import SeqRecord
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

from colored_logging import colored_logging
import logging



def yaml_validator(yaml_file: str, yaml_schema: str) -> None:
    """
    Validates the YAML configuration file by checking for the presence of required fields against provided schema.

        Parameters
        ----------
            :param yaml_file: Path to the YAML configuration file.
            :param yaml_schema: Path to the YAML schema file.

        Returns
        -------
            :returns: bool: True if the YAML file is valid, False otherwise.

        Raises
        ------
            :raises FileNotFoundError: If the specified YAML file does not exist.
            :raises KeyError: If required fields are missing from the YAML.
    """
    try:
        yaml_schema = yamale.make_schema(yaml_schema)
        yaml_data = yamale.make_data(yaml_file)
        yamale.validate(yaml_schema, yaml_data)

        logging.info('YAML configuration file is valid.')

        return True

    except yamale.YamaleError as e:
        for results in e.results:
            for error in results.errors:
                logging.warning(f"Configuration error: {error}")

        return False

def csv_validator(csv_file: str) -> bool:
    """
    Validates the CSV file by checking for the presence of required fields.

        Parameters
        ----------
            :param csv_file: Path to the CSV file.

        Returns
        -------
            :returns: bool: True if the CSV file is valid, False otherwise.

        Raises
        ------
            :raises FileNotFoundError: If the specified CSV file does not exist.
            :raises KeyError: If required fields are missing from the CSV.
    """
    def check_ncbi(series: list) -> bool:
        """
        Checks if the given accession ID exists in the NCBI database.

            Parameters
            ----------
                :param series: The whole accession ID field to check.

            Returns
            -------
                :returns: bool: Series of boolean flags if NCBI entries exists or not.
        """

        def ncbi_exists(acc: str) -> bool:
            time.sleep(0.3)
            logging.debug(f'Checking NCBI for {acc}')
            Entrez.email = defaults.ENTREZ_EMAIL
            try:
                with Entrez.esearch(db='protein', term=acc) as handle:
                    record = Entrez.read(handle)
                    return int(record["Count"]) > 0
            except Exception:
                return False

        return series.apply(lambda x: ncbi_exists(x))

    # Read the CSV file into a DataFrame
    input_df = pd.read_csv(csv_file)

    # Check for required fields
    field_schema = {
        'Label': pa.Column(str, nullable=False),
        'Name': pa.Column(str, nullable=False),
        'Abbreviation': pa.Column(str, nullable=False),
        'Probe': pa.Column(str, nullable=False),
        'Accession': pa.Column(
            str, nullable=False,
            checks=pa.Check(check_ncbi, element_wise=False)
        )
    }

    # Parse schema
    csv_schema = pa.DataFrameSchema(field_schema)

    # Validate
    try:
        csv_schema.validate(input_df)
        logging.info('CSV input file is valid.')
        return True
    except (pa.errors.SchemaError, FileNotFoundError, KeyError) as e:
        logging.warning(f"CSV input error: {e}")
        return False


def fasta_validator(fasta_file: str) -> bool:
    """
    Validates a FASTA file: structure, headers, and (optionally) sequence content.

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
                logging.warning(f"Record {i+1} is missing a header.")
                return False

        logging.info(f'FASTA file {fasta_file} is valid.')
        return True

    except Exception as e:
        return False

def validate_programs() -> bool:
    def check_version(cmd: str, version_cmd: str) -> bool:
        try:
            subprocess.run([cmd, version_cmd], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False

    blast_exists = check_version('tblastn', '-version')
    gt_exists = check_version('gt', '--version')
    datasets_exists = check_version('datasets', '--version')

    if not blast_exists:
        logging.warning('BLAST+ not found. Please install it.')
    else:
        logging.info('BLAST+ is installed.')
    if not gt_exists:
        logging.warning('Genometools not found. Please install it.')
    else:
        logging.info('Genometools is installed.')
    if not datasets_exists:
        logging.warning('Datasets not found. Please install it.')
    else:
        logging.info('Datasets is installed.')

    return all([blast_exists, gt_exists, datasets_exists])


def validate_ncbi_key() -> bool:
    """
    Validates the NCBI API key in the environment. If the key is not found, prompts the user to enter it 
    and sets it in the environment variables.
    """
    if 'NCBI_API_KEY' not in os.environ:
        logging.warning('No NCBI API key found. You can set it with: export NCBI_API_KEY="your_api_key". Do you want to set it now?')
        if api_key := input('Enter your NCBI API key [Leave empty to skip]: ') or None:
            os.environ['NCBI_API_KEY'] = api_key
            logging.info('NCBI API key set in environment variables.')
        else:
            logging.warning('No NCBI API key provided. Expect slower NCBI retrievals.')
    else:
        logging.info("NCBI API key is set in the environment variables.")
        
def main_validator(fasta_files: list) -> bool:
    """
    Main function to validate the YAML, CSV, and FASTA files. It also checks for the NCBI API key in the environment
    and if the required programs are installed.
    
    :return: True if all validations pass, False otherwise.
    """
    logging.debug('Starting input validation process...')

    yaml_validation = yaml_validator(
        yaml_schema=os.path.join(defaults.PATH_DICT['CONFIG_DIR'], 'schema.yaml'),
        yaml_file=defaults.CONFIG_FILE
    )
    
    csv_validation = csv_validator(
        csv_file=defaults.config['input']['probe_csv']
    )

    fasta_bool = []
    for fasta_file in fasta_files:
        fasta_ok = fasta_validator(fasta_file=fasta_file)
        fasta_bool.append(fasta_ok)
    
    fasta_validation = all(fasta_bool)

    program_validation = validate_programs()

    validate_ncbi_key()


    return all([yaml_validation, csv_validation, fasta_validation, program_validation])
        
def green_light(all_valid: bool) -> bool:
    """
    Prints a green light message indicating that the validation is complete.
    """
    if not all_valid:
        logging.warning(f'Settings validation failed. Please check the logs at {defaults.PATH_DICT["LOG_DIR"]} for details.')
    
    if all_valid:
        logging.info('All systems green, ready to rock.')
        time.sleep(0.1)
        proceed = input('Proceed [Y/n]: ') or 'Y'
        if proceed.upper() == 'Y':
            logging.info('RetroSeek started. Depending on your system, this may take some time.')
        elif proceed.upper() == 'N':
            logging.warning('Workflow aborted by user.')
        else:
            logging.warning('Wrong input. Please try again')
            green_light(all_valid)
            
            
if __name__ == '__main__':
    colored_logging(log_file_name='validator.log')
    
    parser.add_argument(
        '--fasta_files',
        nargs='+',
        type=list,
        default=None,
        help='FASTA file paths.'
    )

    args = parser.parse_args()

    # Set FASTA file path list
    fasta_files = args.fasta_files or []

    # Validate all files
    all_valid = main_validator(fasta_files=fasta_files)

    # Print green light message
    green_light(all_valid=all_valid)
    
    