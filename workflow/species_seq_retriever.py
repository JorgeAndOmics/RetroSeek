import os
import io
from io import StringIO
import subprocess
from collections import defaultdict

import defaults
from utils import pickler, unpickler, directory_content_eraser
from object_class import Object
from seq_utils import *
from colored_logging import colored_logging

from Bio.Blast import NCBIXML
from Bio import Entrez

import logging, coloredlogs


def tblastn_retriever(object_dict,
                      command,
                      genome,
                      online_database,
                      input_database_path):
    """
    Orchestrates the tblastn retrieval process. It first performs the tblastn search, then merges the results, and
    finally retrieves the sequences from the online database.

        Args:
            object_dict (dict): A dictionary of objects
            command (str): The BLAST command to run.
            species (str): The species to search for. Retrieved from defaults.
            online_database (str): The online database to retrieve the sequences from.
            input_database_path (str): The path to the local database (species, virus...).

        Returns:
            tblastn_merged2gb_results (dict): A dictionary with the post-BLAST merged and
            retrieved sequences from the online database.
    """
    tblastn_results = blast_threadpool_executor(object_dict=object_dict,
                                                command=command,
                                                genome=genome,
                                                input_database_path=input_database_path)

    tblastn_merged_results = seq_merger(object_dict=tblastn_results)

    tblastn_merged2gb_results = gb_threadpool_executor(object_dict=tblastn_merged_results,
                                                       online_database=online_database)

    return incomplete_dict_cleaner(object_dict=tblastn_merged2gb_results)


if __name__ == '__main__':
    colored_logging(log_file_name='species_seq_retriever.txt')

    probe_dict = unpickler(input_directory_path=os.path.join('..', 'data', 'pickles'),
                           input_file_name='probe_dict.pkl')

    clean_tblastn2gb_results = tblastn_retriever(object_dict=probe_dict,
                                                 command='tblastn',
                                                 genome=defaults.SPECIES,
                                                 online_database='nucleotide',
                                                 input_database_path=defaults.SPECIES_DB)

    directory_content_eraser(directory_path=defaults.TMP_DIR)

    pickler(data=clean_tblastn2gb_results,
            output_directory_path=defaults.PICKLE_DIR,
            output_file_name='tblastn_results.pkl')

    logging.info(f'Successfully performed {command} and retrieval on {len(clean_tblastn_merged2gb_results)} sequences.')