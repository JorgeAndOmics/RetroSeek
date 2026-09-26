"""RetroSeeker_class.py - sequence alignment + GenBank record holder.

This module provides the core functionality for executing BLAST queries against local databases,
parsing the results, and retrieving corresponding GenBank records. It defines a structured RetroSeeker
dataclass for holding sequence alignment information.

Used as part of a larger pipeline for identifying endogenous viral elements (ERVs), this module
supports reproducible, scalable, and programmatic BLAST processing and metadata enrichment.

The ``HSP`` (high-scoring pair) acronym is upper-case throughout this
file because it mirrors BioPython's ``Bio.Blast.Record.HSP`` class
naming. Renaming would diverge the API from the upstream library every
caller is reading. The N802/N803 ruff rules are silenced module-wide
for that reason.

Classes:
- RetroSeeker: stores alignment, HSP, sequence, and metadata for each result.
"""

# ruff: noqa: N802, N803

import tempfile
from dataclasses import dataclass, field
from io import StringIO
from typing import Any

from Bio import SeqIO

import defaults
from log import PipelineError


def _direction(low: Any, high: Any) -> str | None:
    """'+' when ``low < high``, '-' when ``low > high``, None when equal."""
    if low < high:
        return "+"
    if low > high:
        return "-"
    return None


@dataclass
class RetroSeeker:
    """A class to store data from genomic objects.

        Parameters
        ----------
            label: The label of the virus
            virus: The name of the virus
            abbreviation: The abbreviation of the virus
            species: The BLASTed species or agent
            probe: The probe used to identify genomic regions in ERVs
            accession: The accession number of the current contained sequence
            identifier: A random, 6-character string to uniquely identify the object
            strand: The strand of the HSP object
            alignment: Alignment object from BLAST
            HSP: High-scoring pair object from Alignment
            genbank: The Genbank record of the sequence
            fasta: The fasta record of the sequence
            gff: The GFF record of the sequence

    Methods:
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

    label: str | None = field(default=None)
    virus: str | None = field(default=None)
    abbreviation: str | None = field(default=None)
    species: str | None = field(default=None)
    probe: str | None = field(default=None)
    accession: str | None = field(default=None)
    identifier: str | None = field(default=None)
    alignment: Any | None = field(default=None)
    HSP: Any | None = field(default=None)
    genbank: Any | None = field(default=None)
    fasta: Any | None = field(default=None, init=False, repr=False)
    gff: Any | None = field(default=None, init=False, repr=False)
    strand: Any | None = field(default=None, init=False, repr=False)

    def __repr__(self) -> str:
        return (
            f"{self.label},"
            f"{self.virus},"
            f"{self.abbreviation},"
            f"{self.species},"
            f"{self.probe},"
            f"{self.accession},"
            f"{self.identifier}"
        )

    def __hash__(self) -> int:
        return hash(self.identifier)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, RetroSeeker):
            return self.identifier == other.identifier
        return False

    # Static Methods
    @staticmethod
    def extract_fasta_from_genbank(genbank_record: Any) -> str | None:
        """Extracts the FASTA file from the GenBank record.

            Parameters
            ----------
                :param genbank_record: The GenBank record to extract the FASTA from.

        Returns:
            -------
                :returns: The FASTA file content.

        """
        with StringIO() as handle:
            SeqIO.write(genbank_record, handle, "fasta")
            return handle.getvalue()

    @staticmethod
    def extract_gff_from_genbank(genbank_record: Any) -> str | None:
        """Extracts the GFF file from the GenBank record.

            Parameters
            ----------
                :param genbank_record: The GenBank record to extract the GFF from.

        Returns:
            -------
                :returns: The GFF file content.

        """
        with StringIO() as handle:
            for feature in genbank_record.features:
                gff_line = f"{genbank_record.id}\t{feature.type}\t{feature.location}\n"
                handle.write(gff_line)
            return handle.getvalue()

    @staticmethod
    def extract_strand_from_HSP(HSP_obj: Any) -> str | None:
        """Extract the strand ('+' or '-') of an HSP, or None when it cannot be told.

        If the HSP object contains frame information, the strand comes from the
        frame. If the frame is not available, it comes from the orientation of the
        HSP sbjct_start and sbjct_end values.

            Parameters
            ----------
                :param HSP_obj: The HSP object to extract the strand from.

        Returns:
            -------
                :returns: The strand, '+' or '-', or None when it cannot be told.

        """
        # A zero or non-integer frame, or equal sbjct_start and sbjct_end, has no
        # strand to report. The sign is taken before the type check so that an
        # uncomparable frame raises instead of passing as "no strand".
        if HSP_obj.frame:
            last = HSP_obj.frame[-1]
            sign = _direction(0, last)
            return sign if isinstance(last, int) else None
        return _direction(HSP_obj.sbjct_start, HSP_obj.sbjct_end)

    @staticmethod
    def extract_seq2rec(
        seq_obj: Any, obj_type: str, output_type: str = "seqrecord"
    ) -> Any:
        """Parse text (e.g. FASTA) into a SeqRecord object or a temporary file in tmp directory.

            Parameters
            ----------
                :param seq_obj: The text variable to parse (e.g. Instance.fasta, Instance.gff).
                :param obj_type: The Seq object to parse the FASTA text into (e.g. 'fasta', 'gff').
                :param output_type: The type of output to return. Choose 'seqrecord' or 'tempfile'.

        Returns:
            -------
            :returns: A SeqRecord object **or** Path to temporary file containing the parsed FASTA text.

        Raises:
            ------
                :raise ValueError: If an invalid output_type is provided.

        """
        # Create a temporary file to write the FASTA text
        with tempfile.NamedTemporaryFile(
            mode="w",
            delete=False,
            dir=defaults.PATH_DICT["TMP_DIR"],
            suffix=f".{obj_type}",
        ) as tmp_file:
            tmp_file.write(seq_obj)

        # Determine the return type based on output_type parameter
        if output_type == "seqrecord":
            return SeqIO.read(tmp_file.name, obj_type)  # type: ignore[no-untyped-call]
        if output_type == "tempfile":
            return tmp_file.name
        raise ValueError("Invalid output_type. Choose 'seqrecord' or 'tempfile'.")

    # Getters and Setters
    def get_alignment(self) -> Any:
        """Returns the Alignment file associated with the object.

        Returns:
            ----------
                :returns: Alignment or None: The Alignment file content.

        """
        return self.alignment

    def set_alignment(self, alignment_object: object) -> None:
        """Associates an Alignment file with the object.

        Parameters
        ----------
            :param alignment_object: The Alignment object to associate with the object.

        """
        self.alignment = alignment_object

    def get_HSP(self) -> Any:
        """Retrieves the HSP object.

        Returns:
                -------
                    :returns: The HSP object.

        """
        return self.HSP

    def set_HSP(self, HSP_object: object) -> None:
        """Associates an HSP file with the object. Sets strand attribute from the HSP object.

        Parameters
        ----------
            :param HSP_object: The HSP object to associate with the object.

        """
        self.HSP = HSP_object
        self.strand = self.extract_strand_from_HSP(self.HSP)

    def get_genbank(self, output_type: str | None = None) -> Any:
        """Return the GenBank record of the object, or a temporary file holding it.

        If no output_type is provided, it returns the stored GenBank record. If
        output_type is set to 'tempfile', it returns the path to a temporary file
        holding the record's text.

        CAUTION! The GenBank record is already a SeqRecord object. If output_type is set to 'seqrecord', it will raise
        an error. Use only default or 'tempfile' output_type.

            Parameters
            ----------
                :param output_type: Optional(str): The type of output to return. Choose 'seqrecord' or 'tempfile'.

        Returns:
            -------
                :returns: The FASTA file content.
                Raises PipelineError, naming the probe, if its GenBank record was never fetched.

        Raises:
            ------
                :raise Error: If output_type is not 'tempfile'.

        """
        self._require_genbank()
        if output_type:
            return self.extract_seq2rec(
                seq_obj=str(self.genbank), obj_type="genbank", output_type=output_type
            )
        return self.genbank

    def set_genbank(self, genbank_obj: str) -> None:
        """Parse and store a GenBank record, and derive its FASTA and GFF from it.

        The FASTA and GFF are generated through the [extract_fasta_from_genbank]
        and [extract_gff_from_genbank] methods. An unreadable record raises.

            Parameters
            ----------
                :param genbank_obj: The GenBank file to associate with the object.

        """
        # An unreadable record raises: gb_fetcher retries, then reports an ERROR.
        handle = StringIO(genbank_obj)
        self.genbank = SeqIO.read(handle, "genbank")  # type: ignore[no-untyped-call]
        self.fasta = self.extract_fasta_from_genbank(self.genbank)
        self.gff = self.extract_gff_from_genbank(self.genbank)

    def get_fasta(self, output_type: str | None = None) -> Any:
        """Return the FASTA of the object as text, a SeqRecord, or a temporary file.

        If no output_type is provided, it returns the FASTA as a string. If
        output_type is set to 'seqrecord', it returns the FASTA as a SeqRecord
        object. If output_type is set to 'tempfile', it returns the path to a
        temporary file containing the FASTA.

            Parameters
            ----------
                :param output_type: Optional(str): The type of output to return. Choose 'seqrecord' or 'tempfile'.

        Returns:
            -------
                :returns: The FASTA file content.
                Raises PipelineError, naming the probe, if its GenBank record was never fetched.

        """
        self._require_genbank()
        if output_type:
            return self.extract_seq2rec(
                seq_obj=self.fasta, obj_type="fasta", output_type=output_type
            )
        return self.fasta

    def get_gff(self) -> Any:
        """Returns the GFF file associated with the object.

        Returns:
            -------
                :returns: The GFF file content.
                Raises PipelineError, naming the probe, if its GenBank record was never fetched.

        """
        self._require_genbank()
        return self.gff

    def _require_genbank(self) -> None:
        """Stop, naming the probe, when its NCBI record was never fetched."""
        if not self.genbank:
            raise PipelineError(
                f"{self.accession} ({self.probe}) has no GenBank record: its NCBI "
                "fetch failed",
                hint="rerun ./RetroSeek --probe-extractor, and check the accession "
                "in the probe CSV",
            )

    # Display methods and Verifier methods
    def display_info(self) -> str:
        """Return human-readable information about the object as text.

        If an HSP is associated with the object, its information is included too.

        Returns:
            -------
                :returns: str or None: Object information.

        """
        info = (
            f"Label: {self.label}\n"
            f"Virus: {self.virus}\n"
            f"Abbreviation: {self.abbreviation}\n"
            f"Probe: {self.probe}\n"
            f"Accession: {self.accession}\n"
        )

        if self.species:
            info += f"Species: {self.species}\n"

        if self.HSP:
            info += (
                f"Identifier: {self.identifier}\n"
                f"HSP Start: {self.HSP.sbjct_start}\n"
                f"HSP End: {self.HSP.sbjct_end}\n"
                f"HSP Length: {self.HSP.align_length}\n"
                f"HSP Strand: {self.strand}\n"
            )

        if self.is_complete():
            # is_complete() guarantees self.alignment is set; mypy can't narrow
            # through the helper call.
            info += f"Complete Record: {self.alignment.hit_def}\n"  # type: ignore[union-attr]

        return info

    def display_alignment(self) -> str | None:
        """Displays human-readable information about the object's alignment.

        Returns:
            -------
                :returns: str or None: The instance's Alignment information.

        """
        return f"Alignment:\n {self.alignment}\n"

    def display_HSP(self) -> str | None:
        """Displays human-readable information about the object's HSP.

        Returns:
            -------
                :returns: str or None: The instance's HSP information.

        """
        return f"HSP:\n {self.HSP}\n"

    def display_genbank(self) -> str | None:
        """Displays human-readable information about the object's Genbank record.

        Returns:
            -------
                :returns: str or None: The instance's Genbank information.

        """
        return f"Genbank:\n {self.genbank}\n"

    def display_fasta(self) -> str | None:
        """Displays human-readable information about the object's FASTA file.

        Returns:
            -------
                :returns: str or None: The instance's FASTA file content.

        """
        return f"Fasta:\n {self.fasta}\n"

    def display_gff(self) -> str | None:
        """Displays human-readable information about the GFF file associated with the object.

        Returns:
            -------
                :returns: HSP or None: The GFF file content.

        """
        return f"GFF:\n {self.gff}\n"

    def is_complete(self) -> bool:
        """Checks if the object contains Genbank, FASTA and GFF records.

        Returns:
            -------
                :returns: True if the object contains all three records, False otherwise.

        """
        return bool(self.alignment and self.HSP and self.genbank)
