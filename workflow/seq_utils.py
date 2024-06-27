from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO
import pandas as pd
import subprocess
import tempfile
import defaults
import pickle
import random
import string
import time
import os
import re

from object_class import Object
from utils import pickler, unpickler, directory_content_eraser, incomplete_dict_cleaner
from colored_logging import colored_logging

from Bio import SeqIO, Entrez
from Bio.Blast import NCBIXML, NCBIWWW

import logging, coloredlogs


def blaster(instance, command: str, species_database_path, unit: str, outfmt:str= '5'):
    """
    Run a tblastn search for a given object against a given database

        Args:
            instance: The Object instance containing information about the query.
            command: The command to run tblastn.
            species_database_path: The path to the database.
            unit: The particular species against whose database it's being BLASTed
            outfmt: The output format for the BLAST results. Default is 5.

        Returns:
            blast_output: The output of the tblastn search, captured from std_out.

        Raises:
            Exception: If an error occurs while running tblastn.
    """
    with tempfile.NamedTemporaryFile(mode='w', delete=False, dir=defaults.TMP_DIR) as fasta_temp_file:
        if instance.get_fasta():
            fasta_temp_file.write(instance.get_fasta())
            fasta_temp_path = fasta_temp_file.name
    try:
        # Construct the tblastn command
        tblastn_command = [
            command,
            '-db',
            os.path.join(species_database_path, unit, unit),
            '-query',
            fasta_temp_path,  # BLAST+ doesn't take in SeqRecord objects, but files
            '-evalue',
            str(defaults.E_VALUE),
            '-outfmt',
            outfmt,
        ]

        # Run the command
        logging.debug(f'Running tblastn for {instance.probe} against {unit}\n{instance.display_info()}')
        result = subprocess.run(tblastn_command, capture_output=True, text=True)

        # Delete temporal file
        os.remove(fasta_temp_path)

        # Output captured in result.stdout
        blast_output = result.stdout
        if blast_output.strip():
            return blast_output

        logging.error(f'BLAST output is empty for {instance.probe} against {unit}:\n{result.stderr}')
        return None
    except Exception as e:
        logging.error(f'An error occurred while running tblastn: {str(e)}')
        return None


def blaster_parser(result, instance, unit:str):
    """
    Args:
        result: The result of [blaster] function.
        instance: The Object instance containing information about the query.
        unit: The particular species against whose database it's being BLASTed.

    Returns:
        alignment_dict: A dictionary containing the parsed results of the [blaster] function:
        alignment_dict[f'{alignment.id}-{random_string}'] = Object

    Raises:
        Exception: If an error occurs while parsing the BLAST output.
    """
    alignment_dict: dict = {}
    result = StringIO(result)
    try:
        for record in NCBIXML.parse(result):
            for alignment in record.alignments:
                if not record.alignments:
                    logging.warning('No alignments found.')
                    continue
                for hsp in alignment.hsps:

                    regex_pattern = re.compile(defaults.ACCESSION_ID_REGEX)
                    accession_by_regex = regex_pattern.search(alignment.title.split('|')[-1]).group()

                    random_string: str = ''.join(random.choices(string.ascii_uppercase + string.digits, k=6))
                    instance = Object(
                        family=str(instance.family),
                        virus=str(instance.virus),
                        abbreviation=str(instance.abbreviation),
                        species=unit,
                        probe=str(instance.probe),
                        accession=str(accession_by_regex),
                        identifier=instance.identifier or random_string,
                        alignment=alignment,
                        HSP=hsp,
                    )

                    if instance.HSP.align_length >= defaults.PROBE_MIN_LENGTH[instance.probe]:
                        alignment_dict[f'{accession_by_regex}-{random_string}'] = instance
                        logging.info(f'Added {accession_by_regex}-{random_string} to Alignment Dictionary:'
                                     f'\n{instance.display_info()}')


    except Exception as e:
        logging.error(f'Error parsing BLAST output: {e}')
        return None

    return alignment_dict


def blast_task(command, value, input_database_path, unit):
    """
    Run tblastn for the Entrez-retrieved probe sequences against the species database. This function is used as a task
    in the ThreadPoolExecutor

        Args:
            command (str): The type of BLAST to run
            value (Object): The probe instance
            input_database_path (str): The path to the input database (species, virus...)
            unit (str): The species to run tblastn against. Scientific name joined by '_'

        Returns:
            dict: A dictionary containing the tblastn results parsed by [blaster_parser] function

        Raises:
            Exception: If the tblastn process fails

    """
    try:
        if blast_result := blaster(value, command, input_database_path, unit):
            return blaster_parser(blast_result, value, unit)
        logging.warning(f'Could not parse sequences for {value.probe}, {value.virus} against {unit}')
        return None
    except Exception as e:
        logging.error(f'Error running tblastn for {value.probe}, {value.virus} against {unit}: {e}')
        return None


def blast_threadpool_executor(object_dict, command, input_database_path, species):
    """
    Runs BLAST tasks asyncronally using ThreadPoolExecutor

        Args:
            object_dict (dict): A dictionary containing object pairs
            command (str): The type of BLAST to run
            input_database_path (str): The path to the input database (species, virus...)
            species (list): A list of species to run tblastn against. Scientific name joined by '_'

        Returns:
            full_parsed_results (dict): A dictionary containing the parsed BLAST results
    """
    full_parsed_results = {}

    tasks = []
    with ThreadPoolExecutor() as executor:
        for unit in species:
            tasks.extend(
                executor.submit(
                    blast_task, value=value, command=command, input_database_path=input_database_path, unit=unit
                )
                for key, value in object_dict.items()
            )
        for future in as_completed(tasks):
            if result := future.result():
                full_parsed_results |= result

    return full_parsed_results

def gb_fetcher(instance,
               online_database:str,
               _attempt:int=1,
               max_attempts:int=3,
               expand_by:int=0,
               _entrez_email:str=defaults.ENTREZ_EMAIL):
    """
    Fetch the GenBank results for a given sequence and appends it to the objects

        Args:
            instance: The Object instance containing information about the query.
            online_database: The database to fetch the sequence from.
            _attempt: The number of current attempts to fetch the sequence. Default is 1.
            max_attempts: The maximum number of attempts to fetch the sequence. Default is 3.
            expand_by: The number of nucleotides to expand the fetched sequence by at each side. Default is 0.
            _entrez_email: The email to use for the Entrez API. Default is retrieved from defaults.

        Returns:
            object [Optional]: The Object instance with the fetched sequence appended. The function
            updates the input object, but returns the object for convenience.

        Raises:
            Exception: If an error occurs while fetching the sequence.
    """
    Entrez.email = _entrez_email

    kwargs = {
        'db': online_database,
        'id': str(instance.accession),
        'rettype': 'gb',
        'retmode': 'text',
    }
    if instance.HSP:
        kwargs |= {
            'seq_start': max(1, instance.HSP.sbjct_start + expand_by),
            'seq_stop': instance.HSP.sbjct_end + expand_by,
     }
    try:
        with Entrez.efetch(**kwargs) as handle:
            genbank_record = handle.read()
            instance.set_genbank(genbank_record)
            logging.info(f'Fetched GenBank record for:\n{instance.display_info()}')
        return instance
    except Exception as e:
        logging.warning(f'While fetching the genbank record: {str(e)}')

        if _attempt < max_attempts:
            time.sleep(2 ** _attempt)
            logging.warning(f'Retrying... (attempt {_attempt + 1})')
            return gb_fetcher(instance, online_database, _attempt + 1, max_attempts)
        else:
            logging.error(f'Failed to fetch the GenBank record after {max_attempts} attempts.')
            return instance


def gb_threadpool_executor(object_dict, online_database, expand_by:int=0, max_attempts:int=3):
    """
    Fetches GenBank sequences for the objects in the dictionary using ThreadPoolExecutor

        Args:
            object_dict (dict): A dictionary containing object pairs
            online_database (str): The database to retrieve the sequences from
            expand_by (int): The number of nucleotides to expand the fetched sequence by at each side. Default is 0.
            max_attempts (int): The maximum number of attempts to fetch the sequence. Retrieves from defaults.

        Returns:
            full_retrieved_results (dict): A dictionary containing the input objects + the fetched GenBank sequences

        Raises:
            Exception: If an error occurs while fetching the sequences
    """
    full_retrieved_results = {}

    tasks = []
    with ThreadPoolExecutor() as executor:
        tasks.extend(
            executor.submit(
                gb_fetcher,
                instance=value,
                online_database=online_database,
                expand_by=expand_by,
                max_attempts=max_attempts
            )
            for key, value in object_dict.items()
        )
        for future in as_completed(tasks):
            if result := future.result():
                full_retrieved_results[f'{result.accession}-{result.identifier}'] = result

    return full_retrieved_results

def blast2gb2pickle(object_dict: dict,
                    command: str,
                    online_database: str,
                    input_database_path: str,
                    species: list,
                    output_pickle_directory_path,
                    output_pickle_file_name,
                    expand_by:int=defaults.EXPANSION_SIZE,
                    max_attempts:int=defaults.MAX_RETRIEVAL_ATTEMPTS):
    """
    Orchestrates the BLAST process, adds GenBank sequences and pickles

        Args:
            object_dict (dict): A dictionary containing object pairs
            command (str): The type of BLAST to run
            online_database (str): The database to retrieve the sequences from
            input_database_path (str): The path to the input database (species, virus...)
            species (list): A list of species to run tblastn against. Scientific name joined by '_'
            output_pickle_directory_path (str): The path to the output pickle directory
            output_pickle_file_name (str): The name of the output pickle file
            expand_by (int): The number of nucleotides to expand the fetched sequence by at each side. Retrieves from defaults.
            max_attempts (int): The maximum number of attempts to fetch the sequence. Retrieves from defaults.

        Returns:
            None

        Raises:
            None
    """
    blast_results:dict = blast_threadpool_executor(object_dict=object_dict,
                                                command=command,
                                                input_database_path=input_database_path,
                                                species=species)

    if not blast_results:
        logging.critical('BLAST results are empty. Exiting.')
        return

    blast_gb_results:dict = gb_threadpool_executor(object_dict=blast_results,
                                              online_database=online_database,
                                              expand_by=expand_by,
                                              max_attempts=max_attempts)

    if not isinstance(blast_gb_results, dict):
        logging.critical('Fetched GenBank results are not a dictionary. Exiting.')
        return

    logging.debug('Cleaning up incomplete objects...')
    clean_blast_gb_results = incomplete_dict_cleaner(blast_gb_results)

    pickler(clean_blast_gb_results, output_pickle_directory_path, output_pickle_file_name)

    logging.debug('Cleaning up temporal files...')
    directory_content_eraser(defaults.TMP_DIR)

    logging.info(f'Successfully performed {command} and retrieval on {len(blast_gb_results)} sequences.')

