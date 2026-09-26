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
    """One probe, or one BLAST hit of a probe, with its metadata and records.

    Two objects are equal, and hash alike, when their identifiers match. The
    methods fall into four groups: static helpers that derive FASTA, GFF, the
    strand or a SeqRecord from other records; getters and setters (setting the
    HSP also sets the strand, setting the GenBank record also sets the FASTA
    and GFF); display methods that return a record as readable text; and
    is_complete.

    Attributes:
        label: The label of the virus.
        virus: The name of the virus.
        abbreviation: The abbreviation of the virus.
        species: The species the probe was searched against.
        probe: The probe (protein region) used to find ERV regions.
        accession: The accession number of the sequence this object holds.
        identifier: A random 6-character string that identifies the object.
        alignment: The BLAST alignment object.
        HSP: The high-scoring pair object from the alignment.
        genbank: The GenBank record of the sequence, as a SeqRecord.
        fasta: The FASTA text of the sequence, derived from `genbank`.
        gff: The feature table of the sequence, derived from `genbank`.
        strand: The strand of the HSP, '+' or '-', derived from `HSP`.
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
        """Write a GenBank record out as FASTA text.

        Args:
            genbank_record: The GenBank record, as a SeqRecord.

        Returns:
            The record in FASTA format.
        """
        with StringIO() as handle:
            SeqIO.write(genbank_record, handle, "fasta")
            return handle.getvalue()

    @staticmethod
    def extract_gff_from_genbank(genbank_record: Any) -> str | None:
        """List the features of a GenBank record, one tab-separated line each.

        Each line holds the record id, the feature type and the feature
        location. This is a simple feature table, not a full GFF3 file.

        Args:
            genbank_record: The GenBank record, as a SeqRecord.

        Returns:
            The feature lines, or an empty string when the record has none.
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

        Args:
            HSP_obj: The BLAST HSP object.

        Returns:
            The strand, '+' or '-', or None when it cannot be told.
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
        """Turn sequence text (e.g. FASTA) into a SeqRecord or a temporary file.

        The text is always written to a temporary file in the pipeline's
        TMP_DIR first. The file is left in place, whatever the output type.

        Args:
            seq_obj: The text to parse (e.g. the object's `fasta`).
            obj_type: The Biopython format of the text (e.g. 'fasta',
                'genbank'); also used as the file's extension.
            output_type: What to return: 'seqrecord' or 'tempfile'. Defaults
                to 'seqrecord'.

        Returns:
            The text parsed into a single SeqRecord for 'seqrecord', or the
            path of the temporary file, as a string, for 'tempfile'.

        Raises:
            ValueError: If `output_type` is neither 'seqrecord' nor 'tempfile',
                or if the text does not hold exactly one record.
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
        """Return the BLAST alignment object, or None when none is set."""
        return self.alignment

    def set_alignment(self, alignment_object: object) -> None:
        """Store the BLAST alignment object.

        Args:
            alignment_object: The BLAST alignment object.
        """
        self.alignment = alignment_object

    def get_HSP(self) -> Any:
        """Return the BLAST HSP object, or None when none is set."""
        return self.HSP

    def set_HSP(self, HSP_object: object) -> None:
        """Store the BLAST HSP object and set the strand from it.

        Args:
            HSP_object: The BLAST HSP object.
        """
        self.HSP = HSP_object
        self.strand = self.extract_strand_from_HSP(self.HSP)

    def set_genbank(self, genbank_obj: str) -> None:
        """Parse and store a GenBank record, and derive its FASTA and GFF from it.

        The FASTA and GFF are generated through the [extract_fasta_from_genbank]
        and [extract_gff_from_genbank] methods. An unreadable record raises.

        Args:
            genbank_obj: The GenBank record as text, as fetched from Entrez.

        Raises:
            ValueError: If the text does not hold exactly one GenBank record.
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

        Args:
            output_type: None, 'seqrecord' or 'tempfile'.

        Returns:
            The FASTA text, a SeqRecord, or the path of the temporary file.

        Raises:
            PipelineError: If the GenBank record was never fetched; the message
                names the probe.
            ValueError: If `output_type` is not one of the values above.
        """
        self._require_genbank()
        if output_type:
            return self.extract_seq2rec(
                seq_obj=self.fasta, obj_type="fasta", output_type=output_type
            )
        return self.fasta

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

        The species is included when set. The identifier and the HSP
        coordinates, length and strand are included when an HSP is set, and
        the hit definition when the object is complete (see [is_complete]).

        Returns:
            The information, one "Field: value" line each.
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

    def display_fasta(self) -> str | None:
        """Return the FASTA text under a "Fasta:" line."""
        return f"Fasta:\n {self.fasta}\n"

    def display_gff(self) -> str | None:
        """Return the feature table under a "GFF:" line."""
        return f"GFF:\n {self.gff}\n"

    def is_complete(self) -> bool:
        """Tell whether the object holds an alignment, an HSP and a GenBank record.

        Returns:
            True when all three are set, False otherwise.
        """
        return bool(self.alignment and self.HSP and self.genbank)
