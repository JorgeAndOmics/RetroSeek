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

    #TODO: REIMPLEMENT EXISTING DATABASE VERIFICATION AND SKIP

    # Generate the indexes with Suffixerator
    logging.debug('Generating Suffixerator indexes')
    db_utils.ltr_index_generator(input_directory_path=defaults.SPECIES_DB,
                                 file_list=species_files)

    # Generate the LTRHarvest output
    # logging.debug('Generating LTRHarvest output')
    # db_utils.ltr_harvester(index_directory_path=defaults.SPECIES_DB,
    #                        file_list=species_files,
    #                        output_directory_path=defaults.LTRHARVEST_DIR,
    #                        force_rerun=True)
    #
    # # Generate database with LTRHarvest output
    # logging.debug('Generating LTRHarvest database')
    # db_utils.directory_db_generator(file_list=species_files,
    #                                 input_db=defaults.LTRHARVEST_DIR,
    #                                 db_type='nucl',
    #                                 tax_id_input=False,
    #                                 output_directory_path=defaults.LTR_DB)
