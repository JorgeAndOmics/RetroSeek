from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO
import logging

from object_class import Object

from Bio.Blast import NCBIWWW, NCBIXML

import utils
import defaults


def _online_blaster(instance, command: str, online_database: str, _outfmt: str = '5'):
    """
    Runs a BLAST search for a given object against a given database

        Parameters
        ----------
            :param instance: The Object instance containing information about the query.
            :param command: The command to run tblastn.
            :param online_database: The database to run the BLAST against.
            :param _outfmt: The output format for the BLAST results. Default is 5.

        Returns
        -------
            :returns: The output of the BLAST search, captured from std_out.

        Raises
        ------
            :raise Exception: If an error occurs while running BLAST.
    """
    try:
        result_handle = NCBIWWW.qblast(command, online_database, str(instance.get_genbank().seq))
        return NCBIXML.read(result_handle)
    except Exception as e:
        logging.error(f'An error occurred while running {command}: {str(e)}')
        return None


def _online_blaster_parser(result, instance: object) -> dict:
    """
    Parses the output of the [blaster] function.

        Parameters
        ----------
        :param result: The result of [blaster] function.
        :param instance: The Object instance containing information about the query.

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
                        species=instance.species,
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


def _online_blast_task(instance: object, command: str, online_database: str) -> dict:
    """
    Run BLAST command against the online database. This function is used as a task
    in the ThreadPoolExecutor

        Parameters
        ----------
            :param instance: The Object instance containing information about the query.
            :param command: The type of BLAST to run
            :param online_database: The database to run the BLAST against.

        Returns
        -------
            :returns: A dictionary containing the BLAST results parsed by [blaster_parser] function

        Raises
        ------
            :raises Exception: If the BLAST process fails

    """
    try:
        if blast_result := _online_blaster(instance=instance,
                                           command=command,
                                           online_database=online_database):
            return _online_blaster_parser(blast_result, instance)
        logging.warning(f'Could not parse sequences for {instance.probe}, {instance.virus}')
        return None
    except Exception as e:
        logging.error(f'Error running BLAST for {instance.probe}, {instance.virus}: {e}')
        return None


def online_blast_threadpool_executor(object_dict: dict,
                                      command: str,
                                      online_database: str) -> dict:
    """
    Runs BLAST tasks asynchronously against an online database using ThreadPoolExecutor

        Parameters
        ----------
            :param object_dict: A dictionary containing object pairs
            :param command: The type of BLAST to run
            :param online_database: The database to run the BLAST against.

        Returns
        -------
            :returns: A dictionary containing the parsed BLAST results
    """
    full_parsed_results = {}

    tasks = []
    with ThreadPoolExecutor(max_workers=defaults.MAX_THREADPOOL_WORKERS) as executor:

        tasks.extend(
            executor.submit(
                _online_blast_task,
                instance=value,
                command=command,
                online_database=online_database
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
