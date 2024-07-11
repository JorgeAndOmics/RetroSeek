from dataclasses import dataclass, field
from typing import Optional, Any
from io import StringIO
import tempfile
import logging

import defaults

from Bio import SeqIO


@dataclass
class Object:
    """
    A class to store data from genomic objects.

        Parameters
        ----------
            family: The family of the virus
            virus: The name of the virus
            abbreviation: The abbreviation of the virus
            species: The BLASTed species or agent
            probe: The probe used to identify genomic regions in ERVs
            accession: The accession number of the current contained sequence
            identifier: A random, 6-character string to uniquely identify the object
            alignment: Alignment object from BLAST
            HSP: High-scoring pair object from Alignment
            genbank: The Genbank record of the sequence
            fasta: The fasta record of the sequence
            gff: The GFF record of the sequence

        Methods
        -------
            extract_fasta_from_genbank: Extracts the FASTA record from the Genbank record
            extract_gff_from_genbank: Extracts the GFF record from the Genbank record
            extract_strand_from_HSP: Extracts the strand information from the HSP object
            extract_seq2rec: Parse text (eg: FASTA) into a SeqRecord object or a temporary file in tmp directory
            get_alignment: Retrieves the alignment object
            set_alignment: Sets the alignment object
            get_HSP: Retrieves the HSP object
            set_HSP: Sets the HSP object
            get_genbank: Retrieves the Genbank record
            set_genbank: Sets the Genbank record. Also sets the FASTA and GFF records from the Genbank record
            get_fasta: Retrieves the FASTA record
            get_gff: Retrieves the GFF record
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
    identifier: str = field(default=None)
    alignment: Optional[Any] = field(default=None)
    HSP: Optional[Any] = field(default=None)
    genbank: Optional[Any] = field(default=None)
    fasta: Optional[Any] = field(default=None, init=False, repr=False)
    gff: Optional[Any] = field(default=None, init=False, repr=False)
    strand: Optional[Any] = field(default=None, init=False, repr=False)

    def __repr__(self):
        return (f'{self.family},'
                f'{self.virus},'
                f'{self.abbreviation},'
                f'{self.species},'
                f'{self.probe},'
                f'{self.accession},'
                f'{self.identifier}')

    def __hash__(self):  # Needs a hash in order to be used as a dictionary key
        return hash(self.identifier)

    def __eq__(self, other):
        if isinstance(other, Object):
            return self.identifier == other.identifier
        return False

    # Static Methods
    @staticmethod
    def extract_fasta_from_genbank(genbank_record) -> str or None:
        """
        Extracts the FASTA file from the GenBank record.

            Parameters
            ----------
                :param genbank_record: The GenBank record to extract the FASTA from.

            Returns
            -------
                :returns: The FASTA file content.

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
    def extract_gff_from_genbank(genbank_record) -> Optional[str]:
        """
        Extracts the GFF file from the GenBank record.

            Parameters
            ----------
                :param genbank_record: The GenBank record to extract the GFF from.

            Returns
            -------
                :returns: The GFF file content.

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
    def extract_strand_from_HSP(HSP_obj) -> Optional[str]:
        """
        Extracts the strand information from the HSP object. If the HSP object contains frame information, it will
        return the strand based on the frame. If the frame is not available, it will return the strand based on the
        orientation of the HSP sbjct_start and sbjct_end values.

            Parameters
            ----------
                :param HSP_obj: The HSP object to extract the strand from.

            Returns
            -------
                :returns: The strand information (+ or -).

        """
        try:
            if HSP_obj.frame:
                if HSP_obj.frame[-1] > 0 and isinstance(HSP_obj.frame[-1], int):
                    return '+'
                if HSP_obj.frame[-1] < 0 and isinstance(HSP_obj.frame[-1], int):
                    return '-'
            else:
                if HSP_obj.sbjct_start < HSP_obj.sbjct_end:
                    return '+'
                if HSP_obj.sbjct_start > HSP_obj.sbjct_end:
                    return '-'

        except Exception as e:
            logging.warning(f'Could not extract strand: {e}')
            return None

    @staticmethod
    def extract_seq2rec(seq_obj: str, obj_type: str, output_type: str = 'seqrecord'):
        """
        Parse text (e.g. FASTA) into a SeqRecord object or a temporary file in tmp directory.

            Parameters
            ----------
                :param seq_obj: The text variable to parse (e.g. Instance.fasta, Instance.gff).
                :param obj_type: The Seq object to parse the FASTA text into (e.g. 'fasta', 'gff').
                :param output_type: The type of output to return. Choose 'seqrecord' or 'tempfile'.

            Returns
            -------
            :returns: A SeqRecord object **or** Path to temporary file containing the parsed FASTA text.

            Raises
            ------
                :raise ValueError: If an invalid output_type is provided.

        """
        # Create a temporary file to write the FASTA text
        with tempfile.NamedTemporaryFile(mode='w',
                                         delete=False,
                                         dir=defaults.TMP_DIR,
                                         suffix=f'.{obj_type}') as tmp_file:
            tmp_file.write(seq_obj)

        # Determine the return type based on output_type parameter
        if output_type == 'seqrecord':
            return SeqIO.read(tmp_file.name, obj_type)
        elif output_type == 'tempfile':
            return tmp_file.name
        else:
            raise ValueError("Invalid output_type. Choose 'seqrecord' or 'tempfile'.")

    # Getters and Setters
    def get_alignment(self) -> Optional[str]:
        """
        Returns the Alignment file associated with the object.

            Returns
            ----------
                :returns: Alignment or None: The Alignment file content.

        """
        try:
            return self.alignment
        except Exception as e:
            logging.warning(f'Could not retrieve alignment: {e}')

    def set_alignment(self, alignment_object: object) -> None:
        """
        Associates an Alignment file with the object.

            Parameters
            ----------
                :param alignment_object: The Alignment object to associate with the object.

        """
        self.alignment = alignment_object

    def get_HSP(self) -> object:
        """
        Retrieves the HSP object.

                Returns
                -------
                    :returns: The HSP object.

        """
        try:
            return self.HSP
        except Exception as e:
            logging.warning(f'Could not retrieve HSP: {e}')

    def set_HSP(self, HSP_object: object) -> None:

        """
        Associates an HSP file with the object. Sets strand attribute from the HSP object.

            Parameters
            ----------
                :param HSP_object: The HSP object to associate with the object.

        """
        self.HSP = HSP_object
        self.strand = self.extract_strand_from_HSP(self.HSP)

    def get_genbank(self, output_type=None) -> Optional[str]:
        """
        Returns the GenBank file associated with the object. If no output_type is provided, it will return the GenBank
        record as string. If output_type is set to 'tempfile', it will return the path to a temporary file containing
        the FASTA.

        CAUTION! The GenBank record is already a SeqRecord object. If output_type is set to 'seqrecord', it will raise
        an error. Use only default or 'tempfile' output_type.

            Parameters
            ----------
                :param output_type: Optional(str): The type of output to return. Choose 'seqrecord' or 'tempfile'.

            Returns
            -------
                :returns: str or None: The FASTA file content.

            Raises
            ------
                :raise Error: If output_type is not 'tempfile'.

        """
        if self.genbank:
            try:
                if output_type:
                    return self.extract_seq2rec(seq_obj=str(self.genbank),
                                                obj_type='genbank',
                                                output_type=output_type)
                else:
                    return self.genbank

            except Exception as e:
                logging.warning(f'Could not retrieve genbank: {e}')

        else:
            logging.warning('Genbank not set. Could not retrieve genbank.')
            return None

    def set_genbank(self, genbank_obj: str) -> None:
        """
        Associates a GenBank file with the object. Also sets the FASTA and GFF files from the GenBank file
        generated through the [extract_fasta_from_genbank] and [extract_gff_from_genbank] methods.

            Parameters
            ----------
                :param genbank_obj: The GenBank file to associate with the object.

        """
        try:
            handle = StringIO(genbank_obj)
            self.genbank = SeqIO.read(handle, 'genbank')
            self.fasta = self.extract_fasta_from_genbank(self.genbank)
            self.gff = self.extract_gff_from_genbank(self.genbank)
        except Exception as e:
            logging.warning(f'Could not set genbank: {e}')

    def get_fasta(self, output_type=None):
        """
        Returns the FASTA file associated with the object. If no output_type is provided, it will return the FASTA
        as string. If output_type is set to 'seqrecord', it will return the FASTA as a SeqRecord object. If output_type
        is set to 'tempfile', it will return the path to a temporary file containing the FASTA.

            Parameters
            ----------
                :param output_type: Optional(str): The type of output to return. Choose 'seqrecord' or 'tempfile'.

            Returns
            -------
                :returns: str or None: The FASTA file content.

        """
        if self.genbank:
            try:
                if output_type:
                    return self.extract_seq2rec(seq_obj=self.fasta,
                                                obj_type='fasta',
                                                output_type=output_type)
                else:
                    return self.fasta

            except Exception as e:
                logging.warning(f'Could not retrieve fasta: {e}')
        else:
            logging.warning('Genbank not set. Could not retrieve fasta.')
            return None

    def get_gff(self):
        """
        Returns the GFF file associated with the object.

            Returns
            -------
                :returns: str or None: The GFF file content.

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
        it will also display  HSP information.

            Returns
            -------
                :returns: str or None: Object information.

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
            info += f'Complete Record: {self.alignment.hit_def}\n'

        return info

    def display_alignment(self) -> Optional[str]:
        """
        Displays human-readable information about the object's alignment.

            Returns
            -------
                :returns: str or None: The instance's Alignment information.

        """
        return f'Alignment:\n {self.alignment}\n'

    def display_HSP(self) -> Optional[str]:
        """
        Displays human-readable information about the object's HSP.

            Returns
            -------
                :returns: str or None: The instance's HSP information.

        """
        return f'HSP:\n {self.HSP}\n'

    def display_genbank(self) -> Optional[str]:
        """
        Displays human-readable information about the object's Genbank record.

            Returns
            -------
                :returns: str or None: The instance's Genbank information.

        """
        return f'Genbank:\n {self.genbank}\n'

    def display_fasta(self) -> Optional[str]:
        """
        Displays human-readable information about the object's FASTA file.

            Returns
            -------
                :returns: str or None: The instance's FASTA file content.

        """
        return f'Fasta:\n {self.fasta}\n'

    def display_gff(self) -> Optional[str]:
        """
        Displays human-readable information about the GFF file associated with the object.

            Returns
            -------
                :returns: HSP or None: The GFF file content.

        """
        return f'GFF:\n {self.gff}\n'

    def is_complete(self) -> bool:
        """
        Checks if the object contains Genbank, FASTA and GFF records.

            Returns
            -------
                :returns: True if the object contains all three records, False otherwise.

        """
        return bool(self.alignment and self.HSP and self.genbank)
