"""
This script generates the BLAST databases from a collection of (in principle) tax-id named
FASTA files in the specified directories. The databases are generated in their corresponding subdirectories,
the defaults lists are updated with the names of the databases, and the directories are checked for existing
databases to avoid re-generating them.
This script also generates the Suffixerator index, and uses it to generate the LTRHarvest output.

CAUTION: Misnamed files in the input directory will result in duplicated databases.
CAUTION 2: This script relies on the previous generation of the BLAST databases within their directories.

"""
import os
import pickle

import defaults
import utils
import db_utils

from colored_logging import colored_logging
import logging


if __name__ == '__main__':
    # Set up logging
    colored_logging('blast_db_generator.txt')

    # Ensure that the directories exist
    os.makedirs(defaults.SPECIES_DB, exist_ok=True)
    os.makedirs(defaults.VIRUS_DB, exist_ok=True)
    os.makedirs(defaults.LOG_DIR, exist_ok=True)
    os.makedirs(defaults.PICKLE_DIR, exist_ok=True)
    os.makedirs(defaults.TMP_DIR, exist_ok=True)
    os.makedirs(defaults.TABLE_INPUT_DIR, exist_ok=True)
    os.makedirs(defaults.TABLE_OUTPUT_DIR, exist_ok=True)

    # Retrieve files (not directories) in species and virus directories
    logging.debug('Retrieving files from directory')
    species_files: list = utils.directory_file_retriever(defaults.SPECIES_DB)
    virus_files: list = utils.directory_file_retriever(defaults.VIRUS_DB)
    print(species_files)

    #TODO: REIMPLEMENT EXISTING DATABASE VERIFICATION AND SKIP

    # Generate the databases
    logging.debug('Generating BLAST databases')
    db_utils.directory_db_generator(file_list=species_files,
                                    input_db=defaults.SPECIES_DB,
                                    db_type='nucl',
                                    tax_id_input=True,
                                    output_directory_path=defaults.SPECIES_DB)

    db_utils.directory_db_generator(file_list=virus_files,
                                    input_db=defaults.VIRUS_DB,
                                    db_type='nucl',
                                    tax_id_input=True,
                                    output_directory_path=defaults.VIRUS_DB)

