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
    colored_logging('rec_db_generator.txt')

    # Generate necessary folders
    os.makedirs(defaults.A_END_REC_DB, exist_ok=True)
    os.makedirs(defaults.B_END_REC_DB, exist_ok=True)

    # Unpickle the dictionaries
    a_end_dict: dict = utils.unpickler(input_directory_path=defaults.PICKLE_DIR,
                                       input_file_name='tblastn_results.pkl')
    b_end_dict: dict = utils.unpickler(input_directory_path=defaults.PICKLE_DIR,
                                       input_file_name='probe_dict.pkl')

    # Generate the fasta files for the new databases
    db_utils.objdict2fasta(object_dict=a_end_dict,
                           output_directory_path=defaults.A_END_REC_DB,
                           output_file_name='a.fasta')
    db_utils.objdict2fasta(object_dict=b_end_dict,
                           output_directory_path=defaults.B_END_REC_DB,
                           output_file_name='b.fasta')

    # Get the list of files in the directories
    a_end_files = os.listdir(defaults.A_END_REC_DB)
    b_end_files = os.listdir(defaults.B_END_REC_DB)

    # Generate the databases
    db_utils.directory_db_generator(file_list=a_end_files,  # A -> Nucleotides
                                    input_db=defaults.A_END_REC_DB,
                                    db_type='nucl')

    db_utils.directory_db_generator(file_list=b_end_files,  # B -> Proteins
                                    input_db=defaults.B_END_REC_DB,
                                    db_type='prot')
