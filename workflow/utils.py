import pandas as pd
import os
import io
import time
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Optional, Any
import defaults
import pickle

import random
import string

from Bio import SeqIO, Entrez
from Bio.Blast import NCBIXML, NCBIWWW

import logging, coloredlogs

def colored_logging(file_name: str):
    """
    Sets up logging and configures coloredlogs with the custom fields and level styles
    
        Args:
            file_name: The name of the file to save the log in.
        
        Returns:
            None
    """
    # Configure coloredlogs with the custom field and level styles
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(message)s', handlers=[
        logging.FileHandler(os.path.join(defaults.LOG_DIR, file_name), mode='w'),
        logging.StreamHandler()
    ]
                        )

    coloredlogs.install(
        level='DEBUG',
        fmt='%(asctime)s - %(message)s',
        level_styles=defaults.LEVEL_STYLES,
        field_styles=defaults.FIELD_STYLES
    )

def pickler(data, directory_path, file_name):
    """
    Pickles the data to the specified path.

    Args:
        data (Any): The data to be pickled.
        directory_path (str): The directory where the file will be saved.
        file_name (str): The name of the file to save the pickled data in.

    Raises:
        OSError: If the directory cannot be created.
        IOError: If the file cannot be written.
    """
    if not os.path.exists(os.path.join(directory_path)):
        os.makedirs(os.path.join(directory_path), exist_ok=True)
    with open(os.path.join(directory_path, file_name), 'wb') as f:
        pickle.dump(data, f)

def unpickler(directory_path, file_name):
    """
    Unpickles the data from the specified path.

    Args:
        directory_path (str): The directory where the file is saved.
        file_name (str): The name of the file to load the pickled data from.

    Returns:
        Any: The unpickled data.

    Raises:
        FileNotFoundError: If the file cannot be found.
        IOError: If the file cannot be read.
    """
    with open(os.path.join(directory_path, file_name), 'rb') as f:
        return pickle.load(f)


def seq_fetcher(instance, online_database, attempt = 1, max_attempts=3):
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
    logging.info(instance.display_info())
    kwargs = {
        'db': str(online_database),
        'id': instance.accession,
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
    with tempfile.NamedTemporaryFile(mode='w', delete=True, dir=defaults.TMP_DIR) as fasta_temp_file:
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
    result = io.StringIO(result)
    try:
        for record in NCBIXML.parse(result):
            for alignment in record.alignments:
                if not record.alignments:
                    logging.warning('No alignments found.')
                    continue
                for hsp in alignment.hsps:
                    
                    random_string: str = ''.join(random.choices(string.ascii_uppercase + string.digits, k=6))

                    alignment_dict[f'{alignment.hit_id}-{random_string}'] = Object(
                        family=str(object.family),
                        virus=str(object.virus),
                        abbreviation=str(object.abbreviation),
                        species=str(unit),
                        probe=str(object.probe),
                        accession=str(alignment.hit_id),
                        alignment=alignment,
                        HSP=hsp
                    )


    except Exception as e:
        logging.error(f'Error parsing BLAST output: {e}')
        return None

    return alignment_dict

@dataclass
class Object:
    '''
    A class to store data from genomic objects
    '''
    family: str = field(default=None)
    virus: str = field(default=None)
    abbreviation: str = field(default=None)
    species: str = field(default=None)
    probe: str = field(default=None)
    accession: str = field(default=None)
    counter: int = field(default=None)
    alignment: Optional[Any] = field(default=None)
    HSP: Optional[Any] = field(default=None)
    genbank: Optional[Any] = field(default=None)
    fasta: Optional[Any] = field(default=None, init=False, repr=False)
    gff: Optional[Any] = field(default=None, init=False, repr=False)

    def get_alignment(self):
        try:
            return self.alignment
        except Exception as e:
            logging.warning(f'Could not retrieve alignment: {e}')

    def set_alignment(self, alignment):
        self.alignment = alignment

    def get_HSP(self):
        try:
            return self.HSP
        except Exception as e:
            logging.warning(f'Could not retrieve HSP: {e}')

    def set_HSP(self, HSP):
        self.HSP = HSP

    def get_genbank(self):
        try:
            return self.genbank
        except Exception as e:
            logging.warning(f'Could not retrieve genbank: {e}')

    def set_genbank(self, genbank):
        try:
            handle = StringIO(genbank)
            self.genbank = SeqIO.read(handle, 'genbank')
            self.fasta = self.extract_fasta_from_genbank(self.genbank)
            self.gff = self.extract_gff_from_genbank(self.genbank)
        except Exception as e:
            logging.warning(f'Could not set genbank: {e}')

    def get_fasta(self):
        if self.genbank:
            try:
                return self.fasta
            except Exception as e:
                logging.warning(f'Could not retrieve fasta: {e}')
        else:
            logging.warning('Genbank not set. Could not retrieve fasta.')
            return None

    def get_gff(self):
        if self.genbank:
            try:
                return self.gff
            except Exception as e:
                logging.warning(f'Could not retrieve gff: {e}')
        else:
            logging.warning('Genbank not set. Could not retrieve gff.')
            return None

    def extract_fasta_from_genbank(self, genbank_record):
        try:
            with StringIO() as handle:
                SeqIO.write(genbank_record, handle, 'fasta')
                fasta_obj = handle.getvalue()
            return fasta_obj
        except Exception as e:
            logging.warning(f'Could not extract fasta: {e}')
            return None

    def extract_gff_from_genbank(self, genbank_record):
        try:
            with StringIO() as handle:
                for feature in genbank_record.features:
                    gff_line = f'{genbank_record.id}\t{feature.type}\t{feature.location}\n'
                    handle.write(gff_line)
                fasta_obj = handle.getvalue()
            return fasta_obj
        except Exception as e:
            logging.warning(f'Could not extract gff: {e}')
            return None


    def __str__(self):
        return (f'{self.family}, '
               f'{self.virus}, '
               f'{self.abbreviation}, '
               f'{self.species}, '
               f'{self.probe}, '
               f'{self.alignment.id}')

    def display_info(self) -> str:
        return (f'Family: {self.family}\n'
                f'Virus: {self.virus}\n'
                f'Abbreviation: {self.abbreviation}\n'
                f'Species: {self.species}\n'
                f'Probe: {self.probe}\n'
                f'Accession: {self.accession}\n')


    def display_alignment(self) -> str:
        return (f'Alignment:\n {self.alignment}\n')

    def display_HSP(self) -> str:
        return (f'HSP:\n {self.HSP}\n')

    def display_genbank(self) -> str:
        return (f'Genbank:\n {self.genbank}\n')

    def display_fasta(self) -> str:
        return f'Fasta:\n {self.fasta}\n'

    def display_gff(self) -> str:
        return f'GFF:\n {self.gff}\n'

    def is_complete(self) -> bool:
        return bool(self.alignment and self.HSP and self.genbank)





