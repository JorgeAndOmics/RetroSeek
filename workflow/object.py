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



