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

from colored_logging import colored_logging

import logging, coloredlogs

@dataclass
class Object:
    """
    A class to store data from genomic objects

        Attributes:
            family: The family of the virus
            virus: The name of the virus
            abbreviation: The abbreviation of the virus
            species: The BLASTed species or agent
            probe: The probe used to identify genomic regions in ERVs
            accession: The accession number of the current contained sequence
            counter: The number of times the accession has been identified
            identifier: A random, 6-character string to uniquely identify the object
            alignment: Alignment object from BLAST
            HSP: High-scoring pair object from Alignment
            genbank: The Genbank record of the sequence
            fasta: The fasta record of the sequence
            gff: The GFF record of the sequence

        Methods:
            get_alignment: Retrieves the alignment object
            set_alignment: Sets the alignment object
            get_HSP: Retrieves the HSP object
            set_HSP: Sets the HSP object
            get_genbank: Retrieves the Genbank record
            set_genbank: Sets the Genbank record. Also sets the FASTA and GFF records from the Genbank record
            get_fasta: Retrieves the FASTA record
            get_gff: Retrieves the GFF record
            extract_fasta_from_genbank: Extracts the FASTA record from the Genbank record
            extract_gff_from_genbank: Extracts the GFF record from the Genbank record
            display_info: Displays the information contained in the object
            display_alignment: Displays the alignment object
            display_HSP: Displays the HSP object
            display_genbank: Displays the Genbank record
            display_fasta: Displays the FASTA record
            display_gff: Displays the GFF record
            is_complete: Checks if the object contains Genbank, FASTA and GFF records
    """
    family: str = field(default=None)
    virus: str = field(default=None)
    abbreviation: str = field(default=None)
    species: str = field(default=None)
    probe: str = field(default=None)
    accession: str = field(default=None)
    counter: int = field(default=None)
    identifier: str = field(default=None)
    alignment: Optional[Any] = field(default=None)
    HSP: Optional[Any] = field(default=None)
    genbank: Optional[Any] = field(default=None)
    fasta: Optional[Any] = field(default=None, init=False, repr=False)
    gff: Optional[Any] = field(default=None, init=False, repr=False)
    strand: Optional[Any] = field(default=None, init=False, repr=False)

    def __str__(self):
        return (f'{self.family},'
                f'{self.virus},'
                f'{self.abbreviation},'
                f'{self.species},'
                f'{self.probe},'
                f'{self.accession},'
                f'{self.identifier}')

    # Static Methods
    @staticmethod
    def extract_fasta_from_genbank(genbank_record):
        """
        Extracts the FASTA file from the GenBank record

            Returns:
                str or None: The FASTA file content

            Raises:
                None
        """

        try:
            with StringIO() as handle:
                SeqIO.write(genbank_record, handle, 'fasta')
                fasta = handle.getvalue()
            return fasta
        except Exception as e:
            logging.warning(f'Could not extract fasta: {e}')
            return None

    @staticmethod
    def extract_gff_from_genbank(genbank_record):
        """
        Extracts the GFF file from the GenBank record

            Returns:
                str or None: The GFF file content

            Raises:
                None
        """
        try:
            with StringIO() as handle:
                for feature in genbank_record.features:
                    gff_line = f'{genbank_record.id}\t{feature.type}\t{feature.location}\n'
                    handle.write(gff_line)
                gff = handle.getvalue()
            return gff
        except Exception as e:
            logging.warning(f'Could not extract gff: {e}')
            return None

    @staticmethod
    def extract_strand_from_HSP(HSP: Any) -> str:
        """
        Extracts the strand information from the HSP object

            Returns:
                None

            Raises:
                None
        """
        try:
            if HSP.frame:
                if HSP.frame[-1] > 0 and isinstance(HSP.frame[-1], int):
                    return '+'
                if HSP.frame[-1] < 0 and isinstance(HSP.frame[-1], int):
                    return '-'
            else:
                if HSP.sbjct_start < HSP.sbjct_end:
                    return '+'
                if HSP.sbjct_start > HSP.sbjct_end:
                    return '-'

        except Exception as e:
            logging.warning(f'Could not extract strand: {e}')
            return None

    # Getters and Setters
    def get_alignment(self):
        """
        Returns the Alignment file associated with the object

            Returns:
                Alignment or None: The Alignment file content

            Raises:
                None
        """
        try:
            return self.alignment
        except Exception as e:
            logging.warning(f'Could not retrieve alignment: {e}')

    def set_alignment(self, alignment_object):
        """
        Associates an Alignment file with the object

            Returns:
                None

            Raises:
                None
        """
        self.alignment = alignment_object

    def get_HSP(self):
        """
        Retrieves the HSP object

                Arguments: None

                Returns:
                    HSP: The HSP object
        """
        try:
            return self.HSP
        except Exception as e:
            logging.warning(f'Could not retrieve HSP: {e}')

    def set_HSP(self, HSP_object):

        """
        Associates an HSP file with the object. Sets strand attribute from the HSP object.

            Returns:
                None

            Raises:
                None
        """
        self.HSP = HSP_object
        self.strand = self.extract_strand_from_HSP(self.HSP)

    def get_genbank(self):
        """
        Returns the GenBank file associated with the object

            Returns:
                str or None: The GenBank file content

            Raises:
                None
        """

        try:
            return self.genbank
        except Exception as e:
            logging.warning(f'Could not retrieve genbank: {e}')

    def set_genbank(self, genbank):
        """
        Associates a GenBank file with the object. Also sets the FASTA and GFF files from the GenBank file
        generated through the [extract_fasta_from_genbank] and [extract_gff_from_genbank] methods.

            Returns:
                None

            Raises:
                None
        """
        try:
            handle = StringIO(genbank)
            self.genbank = SeqIO.read(handle, 'genbank')
            self.fasta = self.extract_fasta_from_genbank(self.genbank)
            self.gff = self.extract_gff_from_genbank(self.genbank)
        except Exception as e:
            logging.warning(f'Could not set genbank: {e}')

    def get_fasta(self):
        """
        Returns the FASTA file associated with the object

            Returns:
                str or None: The FASTA file content

            Raises:
                None
        """

        if self.genbank:
            try:
                return self.fasta
            except Exception as e:
                logging.warning(f'Could not retrieve fasta: {e}')
        else:
            logging.warning('Genbank not set. Could not retrieve fasta.')
            return None

    def get_gff(self):
        """
        Returns the GFF file associated with the object

            Returns:
                str or None: The GFF file content

            Raises:
                None
        """

        if self.genbank:
            try:
                return self.gff
            except Exception as e:
                logging.warning(f'Could not retrieve gff: {e}')
        else:
            logging.warning('Genbank not set. Could not retrieve gff.')
            return None

    # Display methods and Verifier methods
    def display_info(self) -> str:
        """
        Displays human-readable information about the object. If there are HSPs associated with the object,
        it will also display  HSP information

            Returns:
                str or None: Object information

            Raises:
                None
        """
        info = (f'Family: {self.family}\n'
               f'Virus: {self.virus}\n'
               f'Abbreviation: {self.abbreviation}\n'
               f'Probe: {self.probe}\n'
               f'Accession: {self.accession}\n')

        if self.species:
            info += f'Species: {self.species.replace("_", " ")}\n'

        if self.HSP:
            info += (f'Identifier: {self.identifier}\n'
                     f'HSP Start: {self.HSP.sbjct_start}\n'
                     f'HSP End: {self.HSP.sbjct_end}\n'
                     f'HSP Length: {self.HSP.align_length}\n'
                     f'HSP Strand: {self.strand}\n')

        if self.is_complete():
            info += f'Complete Record'

        return info


    def display_alignment(self) -> str:
        """
        Displays human-readable information about the object's alignment

            Returns:
                str or None: The instance's Alignment information

            Raises:
                None
        """
        return f'Alignment:\n {self.alignment}\n'

    def display_HSP(self) -> str:
        """
        Displays human-readable information about the object's HSP

            Returns:
                str or None: The instance's HSP information

            Raises:
                None
        """
        return f'HSP:\n {self.HSP}\n'

    def display_genbank(self) -> str:
        """
        Displays human-readable information about the object's Genbank record

            Returns:
                str or None: The instance's Genbank information

            Raises:
                None
        """
        return f'Genbank:\n {self.genbank}\n'

    def display_fasta(self) -> str:
        """
        Displays human-readable information about the object's FASTA file

            Returns:
                str or None: The instance's FASTA file content

            Raises:
                None
        """
        return f'Fasta:\n {self.fasta}\n'

    def display_gff(self) -> str:
        """
        Displays human-readable information about the GFF file associated with the object

            Returns:
                HSP or None: The GFF file content

            Raises:
                None
        """
        return f'GFF:\n {self.gff}\n'

    def is_complete(self) -> bool:
        """
        Checks if the object contains Genbank, FASTA and GFF records

            Returns:
                Bool: True if the object contains all three records, False otherwise

            Raises:
                None
        """
        return bool(self.alignment and self.HSP and self.genbank)



