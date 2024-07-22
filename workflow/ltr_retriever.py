"""
This script generates the Suffixerator index, and uses it to generate the LTRHarvest output for the input directory.
"""
import os

import defaults
import utils
import db_utils

from colored_logging import colored_logging
import logging

if __name__ == '__main__':
    # Set up logging
    colored_logging('ltr_retriever.txt')

    # Ensure that the directories exist
    os.makedirs(defaults.SPECIES_DB, exist_ok=True)

    # Retrieve files (not directories) in species and virus directories
    logging.debug('Retrieving files from directory')
    species_files: list = utils.directory_file_retriever(input_directory_path=defaults.SPECIES_DB)

    # Generate a dictionary with the taxid: name pairs
    species_dict: dict = defaults.SPECIES_DICT

    # Clean the file lists if directory with database already exists
    species_files: list = db_utils.existing_ltr_list_cleaner(file_list=species_files,
                                                             directory_to_check=defaults.SPECIES_DB)

    # Generate the indexes with Suffixerator
    logging.debug('Generating Suffixerator indexes')
    db_utils.ltr_index_generator(input_directory_path=defaults.SPECIES_DB,
                                 file_list=species_files,
                                 tax2sc_dict=species_dict)

    # Generate the LTRHarvest output
    logging.debug('Generating LTRHarvest output')
    db_utils.ltr_harvester(input_directory_path=defaults.SPECIES_DB,
                           file_list=species_files,
                           tax2sc_dict=species_dict,
                           output_file_path=defaults.TMP_DIR)
