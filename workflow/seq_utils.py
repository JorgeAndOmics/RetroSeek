from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from io import StringIO
import subprocess
import tempfile
import logging
import random
import string
import time
import os

from object_class import Object

from Bio.Blast import NCBIXML
from Bio import Entrez

import utils
import defaults


def _blaster(instance, command: str, input_database_path, subject: str, _outfmt: str = '5'):
    """
    Runs a BLAST search for a given object against a given database

        Parameters
        ----------
            :param instance: The Object instance containing information about the query.
            :param command: The command to run BLAST.
            :param input_database_path: The path to the database.
            :param subject: The particular genome against whose database it's being BLASTed
            :param _outfmt: The output format for the BLAST results. Default is 5.

        Returns
        -------
            :returns: The output of the BLAST search, captured from std_out.

        Raises
        ------
            :raise Exception: If an error occurs while running BLAST.
    """
    input_path = os.path.join(input_database_path, subject, subject) if subject else input_database_path
    try:
        # Construct the BLAST command
        blast_command = [
            command,
            '-db',
            input_path,
            '-query',
            instance.get_fasta('tempfile'),  # BLAST+ doesn't take in SeqRecord objects, but files
            '-evalue',
            str(defaults.E_VALUE),
            '-outfmt',
            _outfmt,
        ]

        # Run the command
        logging.debug(f'Running {command} for {instance.accession}: {instance.probe} probe against {subject}'
                      f'\n{instance.display_info()}')
        result = subprocess.run(blast_command, capture_output=True, text=True)

        # Output captured in result.stdout
        blast_output = result.stdout
        if blast_output.strip():
            return blast_output

        logging.error(f'BLAST output is empty for {instance.probe} against {subject}:\n{result.stderr}')
        return None
    except Exception as e:
        logging.error(f'An error occurred while running {command}: {str(e)}')
        return None


def _blaster_parser(result, instance: object, subject: str)-> dict:
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
    result = StringIO(result)
    try:
        for record in NCBIXML.parse(result):
            for alignment in record.alignments:
                if not record.alignments:
                    logging.warning('No alignments found.')
                    continue
                for hsp in alignment.hsps:

                    accession_id = alignment.hit_def.split(' ')[0]
                    random_string = utils.random_string_generator(6)

                    new_instance = Object(
                        family=str(instance.family),
                        virus=str(instance.virus),
                        abbreviation=str(instance.abbreviation),
                        species=instance.species or subject,
                        # In order to use the function for virus, species comes already from the input object
                        probe=str(instance.probe),
                        accession=accession_id,
                        identifier=random_string)  # If identifier is inherited from the instance, multiple HSPs will share the same identifier

                    new_instance.set_alignment(alignment),
                    new_instance.set_HSP(hsp),

                    if new_instance.HSP.align_length >= defaults.PROBE_MIN_LENGTH[new_instance.probe]:
                        alignment_dict[f'{accession_id}-{random_string}'] = new_instance
                        logging.info(f'Added {accession_id}-{random_string} to Alignment Dictionary:'
                                     f'\n{new_instance.display_info()}\n')


    except Exception as e:
        logging.error(f'Error parsing BLAST output: {e}')
        return None

    return alignment_dict


def _blast_task(instance: object, command: str, subject: str, input_database_path) -> dict:
    """
    Run BLAST command for the Entrez-retrieved sequences against the species database. This function is used as a task
    in the ThreadPoolExecutor

        Parameters
        ----------
            :param instance: The Object instance containing information about the query.
            :param command: The type of BLAST to run
            :param subject: The species to run BLAST against. Scientific name joined by '_'

        Returns
        -------
            :returns: A dictionary containing the BLAST results parsed by [blaster_parser] function

        Raises
        ------
            :raises Exception: If the BLAST process fails

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
        logging.error(f'Error running BLAST for {instance.probe}, {instance.virus} against {subject.replace("_", " ")}: {e}')
        return None


def blast_threadpool_executor(object_dict: dict,
                              command: str,
                              input_database_path,
                              genome: list = None):
    """
    Runs BLAST tasks asynchronously using ThreadPoolExecutor

        Parameters
        ----------
            :param object_dict: A dictionary containing object pairs
            :param command: The type of BLAST to run
            :param input_database_path: The path to the input database (species, virus...)
            :param genome: Optional: A list of genomes to run BLAST against (Mammals, Virus...), in order to locate the relevant database. Scientific name joined by '_'. If no genome is provided, it just runs the query dictionary against the specified database.


        Returns
        -------
            :returns: A dictionary containing the parsed BLAST results
    """
    full_parsed_results = {}

    tasks = []
    with ThreadPoolExecutor(max_workers=defaults.MAX_THREADPOOL_WORKERS) as executor:

        subjects = genome or [None]

        for subject in subjects:
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
        logging.warning('BLAST results are empty.')
        return

    return full_parsed_results

def blast_monothread_executor(object_dict: dict,
                              command: str,
                              input_database_path,
                              genome: list = None):
    """
    Runs BLAST tasks sequentially without using ThreadPoolExecutor

        Parameters
        ----------
            :param object_dict: A dictionary containing object pairs
            :param command: The type of BLAST to run
            :param input_database_path: The path to the input database (species, virus...)
            :param genome: Optional: A list of genomes to run BLAST against (Mammals, Virus...), in order to locate the relevant database. Scientific name joined by '_'. If no genome is provided, it just runs the query dictionary against the specified database.

        Returns
        -------
            :returns: A dictionary containing the parsed BLAST results
    """
    full_parsed_results = {}

    subjects = genome or [None]

    for subject in subjects:
        for key, value in object_dict.items():
            result = _blast_task(
                instance=value,
                command=command,
                subject=subject,
                input_database_path=input_database_path
            )
            if result:
                full_parsed_results.update(result)

    if not full_parsed_results:
        logging.critical('BLAST results are empty. Exiting.')
        return

    return full_parsed_results



def gb_fetcher(instance: object,
               online_database: str,
               _attempt: int = 1,
               max_attempts: int = 3,
               expand_by: int = 0,
               display_warning: bool = True,
               _entrez_email: str = defaults.ENTREZ_EMAIL,
               _entrez_api_token: str = defaults.NCBI_API_TOKEN):
    """
    Fetch the GenBank results for a given sequence and appends it to the objects.

    CAUTION 1: Expansion and merging of sequences are done simultaneously in this function. New overlaps
    may be created by expanding the sequence.

        Parameters
        ----------
            :param instance: The Object instance containing information about the query.
            :param online_database: The database to fetch the sequence from.
            :param _attempt: The number of current attempts to fetch the sequence. Default is 1.
            :param max_attempts: The maximum number of attempts to fetch the sequence. Default is 3.
            :param expand_by: The number of nucleotides to expand the fetched sequence by at each side. Default is 0.
            :param display_warning: Toggle display of request warning messages. Default is True.
            :param _entrez_email: The email to use for the Entrez API. Default is retrieved from defaults.
            :param _entrez_api_token: The API token to use for the Entrez API. Default is retrieved from defaults.

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
    Entrez.api_key = _entrez_api_token

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
            logging.error(f'Failed to fetch the GenBank record after {max_attempts} attempts.')
            return instance


# noinspection PyTypeChecker
def gb_threadpool_executor(object_dict: dict,
                           online_database: str,
                           display_warning: bool = defaults.DISPLAY_REQUESTS_WARNING,
                           expand_by: int = defaults.EXPANSION_SIZE,
                           max_attempts: int = defaults.MAX_RETRIEVAL_ATTEMPTS):
    """
    Fetches GenBank sequences for the objects in an object dictionary using ThreadPoolExecutor

        Parameters
        ----------
            :param object_dict: A dictionary containing object pairs
            :param online_database: The database to retrieve the sequences from
            :param display_warning: Toggle display of request warning messages - in [_gb_fetcher] -. Default from defaults.
            :param expand_by: The number of nucleotides to expand the fetched sequence by at each side. Default is 0.
            :param max_attempts: The maximum number of attempts to fetch the sequence. Retrieves from defaults.

        Returns
        -------
            :returns: A dictionary containing the input objects + the fetched GenBank sequences

        Raises
        ------
            :raise Exception: If an error occurs while fetching the sequences
    """
    full_retrieved_results: dict = {}

    tasks = []
    with ThreadPoolExecutor() as executor:
        tasks.extend(
            executor.submit(
                utils.execution_limiter,
                func=gb_fetcher,
                instance=value,
                online_database=online_database,
                expand_by=expand_by,
                max_attempts=max_attempts,
                display_warning=display_warning)

            for key, value in object_dict.items()
        )
        for future in as_completed(tasks):
            if result := future.result():
                full_retrieved_results[f'{result.accession}-{result.identifier}'] = result
                logging.info(f'Added {full_retrieved_results[f"{result.accession}-{result.identifier}"].accession}-'
                             f'{full_retrieved_results[f"{result.accession}-{result.identifier}"].identifier} to GenBank Dictionary'
                             f'\n{full_retrieved_results[f"{result.accession}-{result.identifier}"].display_info()}\n')

    if len(full_retrieved_results.items()) == 0:
        logging.critical('No fetched GenBank results. Exiting.')
        return

    return full_retrieved_results


# TODO: Correct gb_monothread_executor function
def gb_monothread_executor(object_dict: dict,
                           online_database: str,
                           display_warning: bool = defaults.DISPLAY_REQUESTS_WARNING,
                           expand_by: int = defaults.EXPANSION_SIZE,
                           max_attempts: int = defaults.MAX_RETRIEVAL_ATTEMPTS) -> dict:
    """
    Fetches GenBank sequences for the objects in an object dictionary using single-thread execution.

        Parameters
        ----------
            :param object_dict: A dictionary containing object pairs.
            :param online_database: The database to retrieve the sequences from.
            :param display_warning: Toggle display of request warning messages - in [_gb_fetcher] -. Default in defaults.
            :param expand_by: The number of nucleotides to expand the fetched sequence by at each side. Default is 0.
            :param max_attempts: The maximum number of attempts to fetch the sequence. Retrieves from defaults.

        Returns
        -------
            :returns: full_retrieved_results: A dictionary containing the input objects + the fetched GenBank sequences

        Raises
        ------
            :raises Exception: If an error occurs while fetching the sequences
    """
    full_retrieved_results = {}

    for key, value in object_dict.items():
        result = gb_fetcher(
            instance=value,
            online_database=online_database,
            expand_by=expand_by,
            max_attempts=max_attempts,
            display_warning=display_warning,
        )
        if result:
            full_retrieved_results[f'{result.accession}-{result.identifier}'] = result
            logging.info(f'Added {result.accession}-{result.identifier} to GenBank Dictionary\n{result.display_info()}')

    if len(full_retrieved_results) == 0:
        logging.critical('No fetched GenBank results. Exiting.')
        return

    return full_retrieved_results



def seq_merger(object_dict: dict):
    """
    Function to merge overlapping sequences. The function groups sequences by species, accession, strand, and virus.
    Then it merges the sequences within each group by updating the coordinates of the merged sequence. The result
    is another dictionary with fewer instances, where the HSP.sbjct_start and HSP.sbjct_end coordinates have been
    updated to reflect the merged sequences, in order to download the full sequence from Entrez.

    CAUTION!: Only HSP.sbjct_start and HSP.sbjct_end are updated. The rest of the attributes are not updated.

        Parameters
        ----------
            :param object_dict: A dictionary with object pairs to merge.

        Returns
        -------
            :returns: A dictionary with the merged sequences.
    """

    # Function to group sequences by species, accession, strand, and virus
    def seq_grouper(object_dict: dict):
        logging.debug('Grouping sequences by species, accession, strand, and virus')
        grouped_sequences = defaultdict(list)

        for key_identifier, instance in object_dict.items():
            group_key = (instance.species, instance.accession, instance.strand, instance.virus)
            grouped_sequences[group_key].append(
                (key_identifier, instance))  # It stores both keys and objects in tuple pairs

        return grouped_sequences

    # Helper function to check if two ranges overlap
    def ranges_overlap(start1: int, end1: int, start2: int, end2: int):
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
            instances.sort(key=lambda x: x[1].HSP.sbjct_end)
        else:
            logging.critical('Critical error: Sequences could not be sorted by their start coordinate. Exiting.')
            return merged_dict

        # Initialize the merged_instance with the first instance in the group
        merged_instance = instances[0][1]

        for _, instance in instances[1:]:
            # Determine the true start and end for the instance
            instance_start = min(instance.HSP.sbjct_start, instance.HSP.sbjct_end)
            instance_end = max(instance.HSP.sbjct_start, instance.HSP.sbjct_end)
            merged_start = min(merged_instance.HSP.sbjct_start, merged_instance.HSP.sbjct_end)
            merged_end = max(merged_instance.HSP.sbjct_start, merged_instance.HSP.sbjct_end)

            if ranges_overlap(instance_start, instance_end, merged_start, merged_end):
                # Merge the sequences by updating the coordinates
                new_start = min(instance_start, merged_start)
                new_end = max(instance_end, merged_end)
                if instance.HSP.sbjct_start < instance.HSP.sbjct_end:
                    merged_instance.HSP.sbjct_start = new_start
                    merged_instance.HSP.sbjct_end = new_end
                else:
                    merged_instance.HSP.sbjct_end = new_start
                    merged_instance.HSP.sbjct_start = new_end
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


def blast_retriever(object_dict: dict,
                    command: str,
                    genome: list,
                    online_database: str,
                    input_database_path,
                    display_warning: bool = defaults.DISPLAY_REQUESTS_WARNING,
                    multi_threading: bool = True) -> dict:
    """
    Orchestrates the blast retrieval process. It first performs the blast search, then merges the results, and
    finally retrieves the sequences from the online database and removes incomplete records.

        Parameters
        ----------
            :param object_dict: A dictionary of objects
            :param command: The BLAST command to run.
            :param genome: The species to search for. Retrieved from defaults.
            :param online_database: The online database to retrieve the sequences from.
            :param input_database_path: The path to the local database (species, virus...).
            :param display_warning: Toggle display of request warning messages. Default from defaults.
            :param multi_threading: Use of multi-threading. Default is True. Execution limiter parameters from defaults.

        Returns
        -------
            :returns: A dictionary with the post-BLAST merged and retrieved sequences from the online database.
    """
    if multi_threading:
        logging.debug('Multithread: Active')
        blast_results: dict = blast_threadpool_executor(object_dict=object_dict,
                                                        command=command,
                                                        genome=genome,
                                                        input_database_path=input_database_path)
    else:
        logging.debug('Multithread: Inactive')
        blast_results: dict = blast_monothread_executor(object_dict=object_dict,
                                                        command=command,
                                                        genome=genome,
                                                        input_database_path=input_database_path)

    blast_merged_results: dict = seq_merger(object_dict=blast_results)

    if multi_threading:
        logging.debug('Multithread: Active')
        blast_merged2gb_results: dict = gb_threadpool_executor(object_dict=blast_merged_results,
                                                               online_database=online_database,
                                                               display_warning=display_warning)
    else:
        logging.debug('Multithread: Inactive')
        blast_merged2gb_results: dict = gb_monothread_executor(object_dict=blast_merged_results,
                                                               online_database=online_database,
                                                               display_warning=display_warning)

    return utils.incomplete_dict_cleaner(object_dict=blast_merged2gb_results)