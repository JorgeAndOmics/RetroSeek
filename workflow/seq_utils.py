"""
seq_utils.py

This module provides utilities to perform BLAST searches on sequence objects,
parse the resulting alignments, fetch associated GenBank records, and organize
retrieved information. The workflow supports local and online databases, allows
threaded execution, and includes failure tolerance for sequence retrieval.

Dependencies:
- Biopython (Bio.Blast, Bio.Entrez)
- tqdm (for progress bars)
- logging
- subprocess (for BLAST calls and format conversions)

Main components:
- `blast_retriever`: Orchestrates the BLAST pipeline.
- `blast_executor`: Runs BLAST tasks sequentially.
- `blaster` / `blaster_parser`: Wraps subprocess-based BLAST call and parses output.
- `gb_fetcher`: Fetches sequence from Entrez given accession ID and alignment.
- `gb_executor`: Fetches sequences for a dictionary of objects.
"""

from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from collections import Counter
from io import StringIO
from tqdm import tqdm
import subprocess
import tempfile
import logging
import random
import string
import time
import sys
import re
import os

from RetroSeeker_class import Object

from Bio.Blast import NCBIXML
from Bio import Entrez

import utils
import defaults


def species_divider(object_dict: dict) -> dict:
    """
    Divides the full_genome_dict into different subdictionaries based on the species contained in the objects

        Parameters
        ----------
            :param object_dict: The dictionary containing the objects to be divided

        Returns
        -------
            :return: A dictionary containing the objects divided by species
    """
    species_dict: dict = defaultdict(dict)
    for key, value in object_dict.items():
        species_dict[value.species][key] = value

    return species_dict


def blaster(instance, command: str, input_database_path, subject: str, num_threads: int,  _outfmt: str = '11'):
    """
    Runs a BLAST search for a given object against a given database

        Parameters
        ----------
            :param instance: The Object instance containing information about the query.
            :param command: The command to run BLAST.
            :param input_database_path: The path to the database.
            :param subject: The particular genome against whose database it's being BLASTed
            :param _outfmt: The output format for the BLAST results. Default is 5.
            :param num_threads: The number of threads to use. Default is 1.

        Returns
        -------
            :returns: The output of the BLAST search, captured from std_out.

        Raises
        ------
            :raise Exception: If an error occurs while running BLAST.
    """
    input_path = os.path.join(input_database_path, subject, subject) if subject else input_database_path
    subject = subject or input_database_path
    try:
        blast_command = [
            command,
            '-db', input_path,
            '-query', instance.get_fasta('tempfile'),
            '-evalue', str(defaults.E_VALUE),
            '-outfmt', _outfmt,
            '-num_threads', str(num_threads),
        ]

        result = subprocess.run(blast_command, capture_output=True, text=True)
        blast_output = result.stdout
        if blast_output.strip():
            return blast_output

        logging.error(f'BLAST (outfmt=11) output is empty for {instance.probe} against {subject}:\n{result.stderr}')
        return None
    except Exception as e:
        logging.error(f'An error occurred while running {command}: {str(e)}')
        return None


def blaster_parser(result, instance: object, subject: str) -> dict:
    """
        Parameters
        ----------
        :param result: The result of [blaster] function.
        :param instance: The Object instance containing information about the query.
        :param subject: The particular genome against whose database it's being BLASTed.

    Returns
    -------
        :returns: A dictionary containing the parsed results of the [blaster] function:
        alignment_dict[f'{alignment.id}-{random_string}'] = Object

    Raises
    ------
        :raise Exception: If an error occurs while parsing the BLAST output.

    CAUTION!: This function is specifically designed to parse the output of the [blaster] function.
    """
    alignment_dict: dict = {}
    regex_pattern = re.compile(defaults.ACCESSION_ID_REGEX)
    try:
        # Write ASN.1 to a temp file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.asn') as tmp_asn:
            tmp_asn.write(result)
            tmp_asn_path = tmp_asn.name

        # Convert ASN.1 to XML
        xml_command = [
            'blast_formatter',
            '-archive', tmp_asn_path,
            '-outfmt', '5'
        ]
        xml_result = subprocess.run(xml_command, capture_output=True, text=True)
        if xml_result.returncode != 0:
            logging.error(f'Error converting ASN to XML: {xml_result.stderr}')
            return None

        xml_string = xml_result.stdout
        xml_handle = StringIO(xml_string)

        # Now parse the XML as before
        for record in NCBIXML.parse(xml_handle):
            for alignment in record.alignments:
                if not record.alignments:
                    logging.warning('No alignments found.')
                    continue
                for hsp in alignment.hsps:
                    accession_id_search = regex_pattern.search(alignment.title)
                    accession_id = accession_id_search[0]
                    random_string = utils.random_string_generator(6)

                    new_instance = Object(
                        family=str(instance.family),
                        virus=str(instance.virus),
                        abbreviation=str(instance.abbreviation),
                        species=instance.species or subject,
                        probe=str(instance.probe).strip(),
                        accession=accession_id,
                        identifier=random_string
                    )

                    new_instance.set_alignment(alignment)
                    new_instance.set_HSP(hsp)

                    alignment_dict[f'{accession_id}-{random_string}'] = new_instance

    except Exception as e:
        # logging.error(f'Error parsing BLAST output: {e}')
        return None

    return alignment_dict


def _blast_task(instance: object, command: str, subject: str, input_database_path, num_threads: int) -> dict:
    """
    Run BLAST command for the Entrez-retrieved sequences against the species database. This function is used as a task
    in the ThreadPoolExecutor

        Parameters
        ----------
            :param instance: The Object instance containing information about the query.
            :param command: The type of BLAST to run
            :param subject: The species to run BLAST against. Scientific name joined by '_'
            :param input_database_path: The path to the database.
            :param num_threads: The number of threads to use. Default is 1.

        Returns
        -------
            :returns: A dictionary containing the BLAST results parsed by [blaster_parser] function

        Raises
        ------
            :raises Exception: If the BLAST process fails

    """
    try:
        if blast_result := blaster(instance=instance,
                                    command=command,
                                    subject=subject,
                                    input_database_path=input_database_path,
                                    num_threads=num_threads):
            return blaster_parser(blast_result, instance, subject)
        logging.warning(f'Could not parse sequences for {instance.probe}, {instance.virus} against {subject}')
        return None
    except Exception as e:
        logging.error(
            f'Error in BLAST task for {instance.probe}, {instance.virus} against {subject}: {e}')
        return None


def blast_executor(object_dict: dict,
                              command: str,
                              input_database_path,
                              num_threads: str,
                              display_full_info: bool,
                              genome: str,
                              ) -> dict:
    """
    Runs BLAST tasks sequentially without using ThreadPoolExecutor

        Parameters
        ----------
            :param object_dict: A dictionary containing object pairs
            :param command: The type of BLAST to run
            :param input_database_path: The path to the input database (species, virus...)
            :param num_threads: The number of threads to use. Default is 1.
            :param display_full_info: Toggle display of full information for each fetched sequence. Default is False.
            :param genome: Optional: A genome to run BLAST against (Mammals, Virus...), in order to locate the relevant database. Scientific name joined by '_'. If no genome is provided, it just runs the query dictionary against the specified database.

        Returns
        -------
            :returns: A dictionary containing the parsed BLAST results
    """
    full_parsed_results = {}

    with tqdm(total=len(object_dict), desc=f'Processing {genome}...') as object_bar:
        for key, value in object_dict.items():
            if result := _blast_task(
                    instance=value,
                    command=command,
                    subject=genome,
                    input_database_path=input_database_path,
                    num_threads=num_threads
            ):
                full_parsed_results |= result
                if display_full_info:
                    key_identifier = f'{value.accession}-{value.identifier}'
                    logging.info(f'Added {key_identifier} to Blast Dictionary\n{value.display_info()}')

            object_bar.update(1)

    if not full_parsed_results:
        logging.critical('BLAST results are empty. Exiting.')
        return

    return full_parsed_results


def blast_retriever(object_dict: dict,
                    command: str,
                    genome: str,
                    input_database_path,
                    num_threads: '1',
                    display_full_info: bool = defaults.DISPLAY_OPERATION_INFO
                    ) -> dict:
    """
    Orchestrates the blast retrieval process. It first performs the blast search, then merges the results, and
    finally retrieves the sequences from the online database and removes incomplete records.

        Parameters
        ----------
            :param object_dict: A dictionary of objects
            :param command: The BLAST command to run.
            :param genome: The species to search for. Retrieved from defaults.
            :param input_database_path: The path to the local database (species, virus...).
            :param num_threads: The number of threads to use. Default is 1.
            :param display_warning: Toggle display of request warning messages. Default from defaults.
            :param display_full_info: Toggle display of full information for each fetched sequence. Default is False.

        Returns
        -------
            :returns: A dictionary with the post-BLAST merged and retrieved sequences from the online database.
    """
    return blast_executor(
        object_dict=object_dict,
        command=command,
        genome=genome,
        num_threads=num_threads,
        input_database_path=input_database_path,
        display_full_info=display_full_info,
        )


def gb_fetcher(instance: object,
               online_database: str,
               _attempt: int = 1,
               max_attempts: int = defaults.MAX_RETRIEVAL_ATTEMPTS,
               display_warning: bool = False,
               _entrez_email: str = defaults.ENTREZ_EMAIL):
    """
    Fetch the GenBank results for a given sequence and appends it to the objects.

        Parameters
        ----------
            :param instance: The Object instance containing information about the query.
            :param online_database: The database to fetch the sequence from.
            :param _attempt: The number of current attempts to fetch the sequence. Default is 1.
            :param max_attempts: The maximum number of attempts to fetch the sequence. Default is 3.
            :param display_warning: Toggle display of request warning messages. Default is True.
            :param _entrez_email: The email to use for the Entrez API. Default is retrieved from defaults.

        Returns
        -------
            :returns: The Object instance with the fetched sequence appended. The function
            returns the instance whether it has been updated or not. If the sequence could not be fetched,
            the instance will be returned as is. If the sequence was fetched, the instance will be updated
            with the GenBank record. See [incomplete_dict_cleaner] function to remove incomplete objects.

        Raises
        ------
            :raise Exception: If an error occurs while fetching the sequence.
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
            genbank_record: object = handle.read()
            instance.set_genbank(genbank_record)
        return instance

    except Exception as e:
        if _attempt < max_attempts:
            time.sleep(2 ** _attempt)
            if display_warning:
                logging.warning(f'While fetching the genbank record: {str(e)}. Retrying... (attempt {_attempt + 1})')
            return gb_fetcher(instance, online_database, _attempt + 1, max_attempts)
        else:
            if display_warning:
                logging.error(f'Failed to fetch the GenBank record after {max_attempts} attempts.')
            return instance


def gb_executor(object_dict: dict,
                           online_database: str,
                           display_warning: bool = defaults.DISPLAY_REQUESTS_WARNING,
                           max_attempts: int = defaults.MAX_RETRIEVAL_ATTEMPTS,
                           display_full_info: bool = False) -> dict:
    """
    Fetches GenBank sequences for the objects in an object dictionary using single-thread execution.

        Parameters
        ----------
            :param object_dict: A dictionary containing object pairs.
            :param online_database: The database to retrieve the sequences from.
            :param display_warning: Toggle display of request warning messages - in [_gb_fetcher] -. Default in defaults.
            :param max_attempts: The maximum number of attempts to fetch the sequence. Retrieves from defaults.
            :param display_full_info: Toggle display of full information for each fetched sequence. Default is False.

        Returns
        -------
            :returns: full_retrieved_results: A dictionary containing the input objects + the fetched GenBank sequences

        Raises
        ------
            :raises Exception: If an error occurs while fetching the sequences
    """
    full_retrieved_results = {}

    with tqdm(total=len(object_dict), desc='Fetching GenBank sequences') as object_bar:
        for key, value in object_dict.items():
            if result := gb_fetcher(
                    instance=value,
                    online_database=online_database,
                    max_attempts=max_attempts,
                    display_warning=display_warning,
            ):
                key_identifier = f'{value.accession}-{value.identifier}'
                full_retrieved_results[key_identifier] = result
                if display_full_info:
                    logging.info(f'Added {key_identifier} to GenBank Dictionary\n{result.display_info()}')
            object_bar.update(1)

    if not full_retrieved_results:
        logging.critical('No fetched GenBank results. Exiting.')
        return

    return full_retrieved_results