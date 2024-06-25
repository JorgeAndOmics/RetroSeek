import os
import io
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import defaults
from utils import *

from Bio.Blast import NCBIXML
from Bio import Entrez

import logging, coloredlogs

def tblastn_task(value, species_database_path, unit):
    try:
        if blast_result := blaster(value, 'tblastn', species_database_path, unit):
            parsed_result = blaster_parser(blast_result, value, unit)
            return parsed_result
        else:
            logging.warning(f'Could not parse sequences for {value.probe}, {value.virus} against {unit}')
            return None
    except Exception as e:
        logging.error(f'Error running tblastn for {value.probe}, {value.virus} against {unit}: {e}')
        return None

def tblastn_er(probe_dict, species_database_path, species):
    '''
    Run tblastn for the Entrez-retrieved probe sequences against the species database

        Args:
            probe_dict (dict): A dictionary containing the probe: object pairs from Entrez
            species_database_path (str): The path to the input species database
            species (list): A list of species to run tblastn against. Scientific name joined by '_'
            online_database (str): The database to retrieve the sequences from
    '''
    full_parsed_results = {}

    tasks = []
    with ThreadPoolExecutor() as executor:
        for unit in species:
            for key, value in probe_dict.items():
                tasks.append(executor.submit(tblastn_task, value, species_database_path, unit))

        for future in as_completed(tasks):
            result = future.result()
            if result:
                full_parsed_results |= result

    return full_parsed_results

def species_seq_retriever(probe_dict,
                          online_database,
                          species_database_path,
                          species,
                          output_pickle_directory_path,
                          output_pickle_file_name):
    '''
    Orchestrates the sequence retrieval process
    '''
    tblastn_results = tblastn_er(probe_dict=probe_dict,
               species_database_path=species_database_path,
               species=species)

    for key, value in tblastn_results.items():
        seq_fetcher(value, online_database)
    pickler(tblastn_results, output_pickle_directory_path, output_pickle_file_name)


if __name__ == '__main__':
    colored_logging('species_seq_retriever.txt')

    probe_dict = unpickler(os.path.join('..', 'data', 'pickles'), 'probe_dict.pkl')

    species_seq_retriever(probe_dict=probe_dict,
                          online_database='Nucleotide',
                          species_database_path=defaults.SPECIES_DB,
                          species=defaults.SPECIES,
                          output_pickle_directory_path=defaults.PICKLE_DIR,
                          output_pickle_file_name='tblastn_results.pkl')
