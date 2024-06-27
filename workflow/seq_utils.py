from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
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


def _blaster(instance, command: str, input_database_path, subject: str, _outfmt:str= '5'):
    """
    Run a tblastn search for a given object against a given database

        Args:
            instance: The Object instance containing information about the query.
            command: The command to run tblastn.
            input_database_path: The path to the database.
            subject: The particular genome against whose database it's being BLASTed
            _outfmt: The output format for the BLAST results. Default is 5.

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
            os.path.join(input_database_path, subject, subject),
            '-query',
            fasta_temp_path,  # BLAST+ doesn't take in SeqRecord objects, but files
            '-evalue',
            str(defaults.E_VALUE),
            '-outfmt',
            _outfmt,
        ]

        # Run the command
        logging.debug(f'Running tblastn for {instance.probe} against {subject}\n{instance.display_info()}')
        result = subprocess.run(tblastn_command, capture_output=True, text=True)

        # Delete temporal file
        os.remove(fasta_temp_path)

        # Output captured in result.stdout
        blast_output = result.stdout
        if blast_output.strip():
            return blast_output

        logging.error(f'BLAST output is empty for {instance.probe} against {subject}:\n{result.stderr}')
        return None
    except Exception as e:
        logging.error(f'An error occurred while running tblastn: {str(e)}')
        return None


def _blaster_parser(result, instance, subject:str):
    """
    Args:
        result: The result of [blaster] function.
        instance: The Object instance containing information about the query.
        subject: The particular genome against whose database it's being BLASTed.

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
                    new_instance = Object(
                        family=str(instance.family),
                        virus=str(instance.virus),
                        abbreviation=str(instance.abbreviation),
                        species=instance.species or subject, # In order to use the function for virus, species comes already from the input object
                        probe=str(instance.probe),
                        accession=str(accession_by_regex),
                        identifier=instance.identifier or random_string)

                    new_instance.set_alignment(alignment),
                    new_instance.set_HSP(hsp),


                    if new_instance.HSP.align_length >= defaults.PROBE_MIN_LENGTH[new_instance.probe]:
                        alignment_dict[f'{accession_by_regex}-{random_string}'] = new_instance
                        logging.info(f'Added {accession_by_regex}-{random_string} to Alignment Dictionary:'
                                     f'\n{new_instance.display_info()}')


    except Exception as e:
        logging.error(f'Error parsing BLAST output: {e}')
        return None

    return alignment_dict


def _blast_task(instance, command, subject, input_database_path):
    """
    Run tblastn for the Entrez-retrieved probe sequences against the species database. This function is used as a task
    in the ThreadPoolExecutor

        Args:
            instance: The Object instance containing information about the query.
            command (str): The type of BLAST to run
            subject (str): The species to run tblastn against. Scientific name joined by '_'

        Returns:
            dict: A dictionary containing the tblastn results parsed by [blaster_parser] function

        Raises:
            Exception: If the tblastn process fails

    """
    try:
        if blast_result := _blaster(instance=instance,
                                    command=command,
                                    subject=subject,
                                    input_database_path=input_database_path):
            return _blaster_parser(blast_result, instance, subject)
        logging.warning(f'Could not parse sequences for {instance.probe}, {instance.virus} against {subject}')
        return None
    except Exception as e:
        logging.error(f'Error running tblastn for {instance.probe}, {instance.virus} against {subject}: {e}')
        return None


def blast_threadpool_executor(object_dict,
                              command,
                              genome,
                              input_database_path):
    """
    Runs BLAST tasks asyncronally using ThreadPoolExecutor

        Args:
            object_dict (dict): A dictionary containing object pairs
            command (str): The type of BLAST to run
            genome (list): A list of genome to run tblastn against (Mammals, Virus...). Scientific name joined by '_'.
            input_database_path (str): The path to the input database (species, virus...)


        Returns:
            full_parsed_results (dict): A dictionary containing the parsed BLAST results
    """
    full_parsed_results = {}

    tasks = []
    with ThreadPoolExecutor() as executor:
        for subject in genome:
            tasks.extend(
                executor.submit(
                    _blast_task,
                    instance=value,
                    command=command,
                    subject=subject,
                    input_database_path=input_database_path
                )
                for key, value in object_dict.items()
            )
        for future in as_completed(tasks):
            if result := future.result():
                full_parsed_results |= result

    if not full_parsed_results:
        logging.critical('BLAST results are empty. Exiting.')
        return

    return full_parsed_results

def _gb_fetcher(instance,
               online_database:str,
               _attempt:int=1,
               max_attempts:int=3,
               expand_by:int=0,
               _entrez_email:str=defaults.ENTREZ_EMAIL):
    """
    Fetch the GenBank results for a given sequence and appends it to the objects.

    CAUTION 1: Expansion and merging of sequences are done simultaneously in this function. New overlaps
    may be created by expanding the sequence.

        Args:
            instance: The Object instance containing information about the query.
            online_database: The database to fetch the sequence from.
            _attempt: The number of current attempts to fetch the sequence. Default is 1.
            max_attempts: The maximum number of attempts to fetch the sequence. Default is 3.
            expand_by: The number of nucleotides to expand the fetched sequence by at each side. Default is 0.
            _entrez_email: The email to use for the Entrez API. Default is retrieved from defaults.

        Returns:
            object: The Object instance with the fetched sequence appended. The function
            returns the instance whether it has been updated or not. If the sequence could not be fetched,
            the instance will be returned as is. If the sequence was fetched, the instance will be updated
            with the GenBank record. See [incomplete_dict_cleaner] function to remove incomplete objects.

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
            return _gb_fetcher(instance, online_database, _attempt + 1, max_attempts)
        else:
            logging.error(f'Failed to fetch the GenBank record after {max_attempts} attempts.')
            return instance


def gb_threadpool_executor(object_dict,
                           online_database,
                           expand_by:int=defaults.EXPANSION_SIZE,
                           max_attempts:int=defaults.MAX_RETRIEVAL_ATTEMPTS):
    """
    Fetches GenBank sequences for the objects in an object dictionary using ThreadPoolExecutor

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
                _gb_fetcher,
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
                logging.info(f'Added {result.accession}-{result.identifier} to GenBank Dictionary\n')

    if not full_retrieved_results:
        logging.critical('No fetched GenBank results. Exiting.')
        return

    return full_retrieved_results


# Function to merge overlapping sequences within each group
def seq_merger(object_dict):
    """
    Function to merge overlapping sequences. The function groups sequences by species, accession, strand, and virus.
    Then it merges the sequences within each group by updating the coordinates of the merged sequence. The result
    is another dictionary with fewer instances, where the HSP.sbjct_start and HSP.sbjct_end coordinates have been
    updated to reflect the merged sequences, in order to download the full sequence from Entrez.

    CAUTION 1!: This function assumes that the sequences are sorted by their start coordinate.
    CAUTION 2!: Only HSP.sbjct_start and HSP.sbjct_end are updated. The rest of the attributes are not updated.

        Args:
            object_dict (dict): A dictionary with object pairs to merge.

        Returns:
            merged_dict (dict): A dictionary with the merged sequences.
    """
    # Function to group sequences by species, accession, strand, and virus
    def seq_grouper(object_dict):
        logging.debug('Grouping sequences by species, accession, strand, and virus')
        grouped_sequences = defaultdict(list)

        for key_identifier, instance in object_dict.items():
            group_key = (instance.species, instance.accession, instance.strand, instance.virus)
            grouped_sequences[group_key].append((key_identifier, instance)) # It stores both keys and objects in tuple pairs

        return grouped_sequences

    # Helper function to check if two ranges overlap
    def ranges_overlap(start1, end1, start2, end2):
        return max(start1, start2) <= min(end1, end2)

    grouped_sequences = seq_grouper(object_dict=object_dict)

    merged_dict = {}

    for group_key, instances in grouped_sequences.items():
        if not instances:
            continue

        # Sort instances based on their coordinates
        if all(instance.HSP.sbjct_start < instance.HSP.sbjct_end for _, instance in instances):
            instances.sort(key=lambda x: x[1].HSP.sbjct_start)
        elif all(instance.HSP.sbjct_end < instance.HSP.sbjct_start for _, instance in instances):
            instances.sort(key=lambda x: x[1].HSP.sbjct_start, reverse=True)
        else:
            logging.critical('Critical error: Sequences could not be sorted by their start coordinate. Exiting.')
            return merged_dict

        # Initialize the merged_instance with the first instance in the group
        merged_instance = instances[0][1]

        for _, instance in instances[1:]:
            if ranges_overlap(instance.HSP.sbjct_start, instance.HSP.sbjct_end,
                              merged_instance.HSP.sbjct_start, merged_instance.HSP.sbjct_end):
                # Merge the sequences by updating the coordinates
                merged_instance.HSP.sbjct_start = min(instance.HSP.sbjct_start, merged_instance.HSP.sbjct_start)
                merged_instance.HSP.sbjct_end = max(instance.HSP.sbjct_end, merged_instance.HSP.sbjct_end)
            else:
                # If they do not overlap, add the merged_instance to the dictionary
                new_key = f"{merged_instance.accession}-{merged_instance.identifier}"
                merged_dict[new_key] = merged_instance
                logging.info(f'Added {new_key} to Merging Dictionary')
                # Set the current instance as the new merged_instance
                merged_instance = instance

        # Add the last merged_instance to the dictionary
        new_key = f"{merged_instance.accession}-{merged_instance.identifier}"
        merged_dict[new_key] = merged_instance
        logging.info(f'Added {new_key} to Merging Dictionary')

    return merged_dict