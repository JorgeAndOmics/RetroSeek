from dataclasses import dataclass, field
from typing import Optional, Any
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
        os.unlink(fasta_temp_path)

        # Output captured in result.stdout
        blast_output = result.stdout
        if blast_output.strip():
            return blast_output

        logging.error(f'BLAST output is empty for {instance.probe} against {unit}:\n{result.stderr}')
        return None
    except Exception as e:
        logging.error(f'An error occurred while running tblastn: {str(e)}')
        return None


def blaster_parser(result, object, unit):
    """
    Args:
        result: The result of [blaster] function.
        object: The Object instance containing information about the query.
        unit: The particular species against whose database it's being BLASTed

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
                        family=str(object.family),
                        virus=str(object.virus),
                        abbreviation=str(object.abbreviation),
                        species=str(unit),
                        probe=str(object.probe),
                        accession=str(accession_by_regex),
                        alignment=alignment,
                        HSP=hsp
                    )

                    if instance.HSP.align_length >= defaults.PROBE_MIN_LENGTH[object.probe]:
                        alignment_dict[f'{accession_by_regex}-{random_string}'] = instance
                        logging.info(f'Added {accession_by_regex}-{random_string} -> {instance.display_info()}')


    except Exception as e:
        logging.error(f'Error parsing BLAST output: {e}')
        return None

    return alignment_dict


def seq_fetcher(instance, online_database, attempt=1, max_attempts=3):
    '''
    Fetch the tblastn results for a given sequence and appends it to the objects

        Args:
            instance: The Object instance containing information about the query.
            online_database: The database to fetch the sequence from.
            attempt: The number of current attempts to fetch the sequence. Default is 1.
            max_attempts: The maximum number of attempts to fetch the sequence. Default is 3.

        Returns:
            object [Optional]: The Object instance with the fetched sequence appended. The function
            updates the input object, but returns the object for convenience.

        Raises:
            Exception: If an error occurs while fetching the sequence.
    '''
    Entrez.email = defaults.ENTREZ_EMAIL

    kwargs = {
        'db': str(online_database),
        'id': str(instance.accession),
        'rettype': 'gb',
        'retmode': 'text'
    }
    if instance.HSP:
        kwargs |= {
            'seq_start': instance.HSP.sbjct_start,
            'seq_stop': instance.HSP.sbjct_end,
     }
    try:
        with Entrez.efetch(**kwargs) as handle:
            genbank_record = handle.read()
            instance.set_genbank(genbank_record)
            logging.info(f'Fetched GenBank record for:\n{instance.display_info()}')
        return instance
    except Exception as e:
        logging.error(f'An error occurred while fetching the genbank record: {str(e)}')

        if attempt < max_attempts:
            time.sleep(2 ** attempt)
            logging.info(f'Retrying... (attempt {attempt + 1})')
            return seq_fetcher(instance, online_database, attempt + 1, max_attempts)
        else:
            logging.error(f'Failed to fetch the GenBank record after {max_attempts} attempts.')
            return instance

