"""
This script generates the BLAST databases from a collection of (in principle) tax-id named
FASTA files in the specified directories. The databases are generated in their corresponding subdirectories,
the defaults lists are updated with the names of the databases, and the directories are checked for existing
databases to avoid re-generating them.

CAUTION: Misnamed files in the input directory will result in duplicated databases.
"""
import os

import defaults
import utils
import db_utils

from colored_logging import colored_logging


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

    # Retrieve files in species and virus directories
    species_files = utils.directory_file_retriever(defaults.SPECIES_DB)
    virus_files = utils.directory_file_retriever(defaults.VIRUS_DB)

    # Populate the species and virus lists in defaults
    defaults.SPECIES = list(set(db_utils.taxlist2name(os.listdir(defaults.SPECIES_DB))))
    defaults.VIRUS = list(set(db_utils.taxlist2name(os.listdir(defaults.VIRUS_DB))))

    # Clean the file lists if directory with database already exists
    species_files = db_utils.existing_db_list_cleaner(species_files, defaults.SPECIES_DB)
    virus_files = db_utils.existing_db_list_cleaner(virus_files, defaults.VIRUS_DB)

    # Generate the databases
    db_utils.directory_db_generator(file_list=species_files,
                                    input_db=defaults.SPECIES_DB,
                                    db_type='nucl')

    db_utils.directory_db_generator(file_list=virus_files,
                                    input_db=defaults.VIRUS_DB,
                                    db_type='nucl')
