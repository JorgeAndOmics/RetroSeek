"""BLAST searches of probe objects, hit parsing and GenBank retrieval.

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

import logging
import tempfile
import time
from io import StringIO
from pathlib import Path
from typing import Any

from Bio import Entrez
from Bio.Blast import NCBIXML
from tqdm import tqdm

import defaults
import utils
from external import run_tool
from log import PipelineError
from RetroSeeker_class import RetroSeeker

logger = logging.getLogger(__name__)


def blaster(
    instance: RetroSeeker,
    command: str,
    input_database_path: str | Path,
    subject: str,
    num_threads: int,
    _outfmt: str = "11",
) -> str:
    """Run a BLAST search of one probe object against a genome database.

    The output is the ASN.1 archive (outfmt 11) that blaster_parser reads. An
    empty output stops the job with a PipelineError, since the archive is
    written even when nothing matches.

    Args:
        instance: The probe to search with; its sequence is the BLAST query.
        command: The BLAST program to run, e.g. "tblastn".
        input_database_path: The folder that holds one BLAST database per
            genome, or the database itself when `subject` is empty.
        subject: The genome whose database is searched, which is also the name
            of its database under `input_database_path`.
        num_threads: The number of threads BLAST may use.
        _outfmt: The BLAST output format. Defaults to "11" (ASN.1 archive).

    Returns:
        The BLAST output, captured from standard output.

    Raises:
        PipelineError: If BLAST is missing, fails, or gives no output.
    """
    input_path = (
        str(Path(input_database_path) / subject / subject)
        if subject
        else input_database_path
    )
    subject = subject or str(input_database_path)
    blast_command = [
        command,
        "-db",
        input_path,
        "-query",
        instance.get_fasta("tempfile"),
        "-evalue",
        str(defaults.E_VALUE),
        "-outfmt",
        _outfmt,
        "-num_threads",
        str(num_threads),
    ]
    blast_output = run_tool(blast_command).stdout
    # The ASN.1 archive (outfmt 11) is written even when nothing matches, so an
    # empty one means BLAST broke, not that the probe has no hits.
    if not blast_output.strip():
        raise PipelineError(
            f"{command} gave no output for probe {instance.probe} against {subject}",
            hint="check the BLAST database with --blast-dbs; the job log has the details",
        )
    return blast_output


def _accession(hit_def: str | None) -> str:
    """The accession of a BLAST hit: the first whitespace token of its hit_def.

    FASTA-header convention: "CM138268.1 Molossus molossus chr 3, whole genome
    shotgun sequence" -> "CM138268.1". Storing only the accession keeps the
    seqid identical to how LTRdigest / rtracklayer / GRanges represent it, so
    downstream findOverlaps matches without per-stage stripping.
    """
    raw_hit_def = hit_def or ""
    return raw_hit_def.split()[0] if raw_hit_def else raw_hit_def


def _unused_identifier(used: set[str]) -> str:
    """A random 6-character identifier not yet given to a hit of this genome.

    Every probe of a genome shares one key space ({accession}-{identifier}), and
    RetroSeeker objects compare equal by identifier, so a repeated draw would
    either replace an earlier hit or make two hits indistinguishable. The new
    identifier is recorded in ``used``.
    """
    while True:
        identifier = utils.random_string_generator(6)
        if identifier not in used:
            used.add(identifier)
            return identifier


def _hit_object(
    instance: RetroSeeker,
    subject: str,
    alignment: Any,
    hsp: Any,
    used: set[str],
) -> tuple[str, RetroSeeker]:
    """One HSP as a RetroSeeker carrying the query's metadata, and its dict key.

    The key is ``{accession}-{identifier}``; the identifier is new to ``used``.
    """
    accession_id = _accession(alignment.hit_def)
    random_string = _unused_identifier(used)
    new_instance = RetroSeeker(
        label=str(instance.label),
        virus=str(instance.virus),
        abbreviation=str(instance.abbreviation),
        species=instance.species or subject,
        probe=str(instance.probe).strip(),
        accession=accession_id,
        identifier=random_string,
    )
    new_instance.set_alignment(alignment)
    new_instance.set_HSP(hsp)
    return f"{accession_id}-{random_string}", new_instance


def blaster_parser(
    result: str,
    instance: RetroSeeker,
    subject: str,
    used_identifiers: set[str] | None = None,
) -> dict[str, RetroSeeker] | None:
    """Parse a BLAST archive into one RetroSeeker object per HSP.

    The archive is converted to XML with blast_formatter; each HSP becomes a
    RetroSeeker carrying the query's metadata, keyed by
    ``{accession}-{identifier}``.

    CAUTION: This function is designed to parse only the output of [blaster].

    Args:
        result: The ASN.1 archive returned by [blaster].
        instance: The probe that was searched; its metadata is copied to each
            hit.
        subject: The genome that was searched; used as the hit's species when
            the probe has none.
        used_identifiers: Identifiers already given to this genome's hits; new
            ones are added. Pass the same set for every probe of a genome.

    Returns:
        A dictionary with one RetroSeeker per HSP, keyed by
        ``{accession}-{identifier}``. It is empty when there are no hits.

    Raises:
        PipelineError: If blast_formatter cannot read the archive.
    """
    alignment_dict: dict[str, RetroSeeker] = {}
    used = set() if used_identifiers is None else used_identifiers
    # blast_formatter reads the ASN.1 archive from a file; the file is removed
    # whatever happens. Any failure below stops the job: a half-parsed genome
    # would otherwise look like one with fewer hits.
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".asn") as tmp_asn:
        tmp_asn.write(result)
        tmp_asn_path = tmp_asn.name
    try:
        xml_command = ["blast_formatter", "-archive", tmp_asn_path, "-outfmt", "5"]
        xml_handle = StringIO(run_tool(xml_command).stdout)
        for record in NCBIXML.parse(xml_handle):  # type: ignore[no-untyped-call]
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    key, hit = _hit_object(instance, subject, alignment, hsp, used)
                    alignment_dict[key] = hit
    finally:
        Path(tmp_asn_path).unlink(missing_ok=True)

    return alignment_dict


def _blast_task(
    instance: RetroSeeker,
    command: str,
    subject: str,
    input_database_path: str | Path,
    num_threads: int,
    used_identifiers: set[str] | None = None,
) -> dict[str, RetroSeeker] | None:
    """Run one probe against the species database and parse its hits.

    One task of blast_executor: [blaster] followed by [blaster_parser].

    Args:
        instance: The probe to search with.
        command: The BLAST program to run, e.g. "tblastn".
        subject: The genome to search: its scientific name joined by "_".
        input_database_path: The folder that holds one BLAST database per
            genome.
        num_threads: The number of threads BLAST may use.
        used_identifiers: Identifiers already given to this genome's hits; new
            ones are added.

    Returns:
        The probe's hits as parsed by [blaster_parser].

    Raises:
        PipelineError: If BLAST or blast_formatter fails.
    """
    blast_result = blaster(
        instance=instance,
        command=command,
        subject=subject,
        input_database_path=input_database_path,
        num_threads=num_threads,
    )
    return blaster_parser(blast_result, instance, subject, used_identifiers)


def blast_executor(
    object_dict: dict[str, RetroSeeker],
    command: str,
    input_database_path: str | Path,
    num_threads: int,
    genome: str,
) -> dict[str, RetroSeeker] | None:
    """Run every probe against one genome in turn and merge their hits.

    All probes of the genome share one set of hit identifiers. Returns None
    (after a CRITICAL log line) when no probe has a hit.

    Args:
        object_dict: The probes to search with, as RetroSeeker objects.
        command: The BLAST program to run, e.g. "tblastn".
        input_database_path: The folder that holds one BLAST database per
            genome, or the database itself when `genome` is empty.
        num_threads: The number of threads BLAST may use.
        genome: The genome to search: its scientific name joined by "_", which
            locates its database. When empty, the probes are searched against
            `input_database_path` directly.

    Returns:
        All hits of all probes, keyed by ``{accession}-{identifier}``, or None
        when there are none.

    Raises:
        PipelineError: If BLAST or blast_formatter fails.
    """
    full_parsed_results: dict[str, RetroSeeker] = {}
    used_identifiers: set[str] = set()  # one key space for all probes of a genome

    # disable=None: no bar when stderr is not a terminal (under the launcher),
    # where its carriage returns would flood the run log.
    with tqdm(
        total=len(object_dict), desc=f"Processing {genome}...", disable=None
    ) as object_bar:
        for value in object_dict.values():
            if result := _blast_task(
                instance=value,
                command=command,
                subject=genome,
                input_database_path=input_database_path,
                num_threads=num_threads,
                used_identifiers=used_identifiers,
            ):
                full_parsed_results |= result
                key_identifier = f"{value.accession}-{value.identifier}"
                logger.debug(
                    f"Added {key_identifier} to Blast Dictionary\n{value.display_info()}"
                )

            object_bar.update(1)

    if not full_parsed_results:
        logger.critical("BLAST results are empty. Exiting.")
        return None

    return full_parsed_results


def blast_retriever(
    object_dict: dict[str, RetroSeeker],
    command: str,
    genome: str,
    input_database_path: str | Path,
    num_threads: int,
) -> dict[str, RetroSeeker] | None:
    """Run the BLAST search of every probe against one genome and merge the hits.

    It delegates to [blast_executor]; retrieving GenBank records is done
    separately by [gb_executor].

    Args:
        object_dict: The probes to search with, as RetroSeeker objects.
        command: The BLAST program to run, e.g. "tblastn".
        genome: The genome to search: its scientific name joined by "_".
        input_database_path: The folder that holds one BLAST database per
            genome.
        num_threads: The number of threads BLAST may use.

    Returns:
        All hits of all probes, keyed by ``{accession}-{identifier}``, or None
        when there are none.

    Raises:
        PipelineError: If BLAST or blast_formatter fails.
    """
    return blast_executor(
        object_dict=object_dict,
        command=command,
        genome=genome,
        num_threads=num_threads,
        input_database_path=input_database_path,
    )


def gb_fetcher(
    instance: RetroSeeker,
    online_database: str,
    _attempt: int = 1,
    max_attempts: int = defaults.MAX_RETRIEVAL_ATTEMPTS,
    _entrez_email: str = defaults.ENTREZ_EMAIL,
) -> RetroSeeker:
    """Fetch the GenBank record for one object from Entrez and attach it.

    When the object has an HSP, only the aligned stretch of the subject
    sequence is fetched. A failed fetch is retried after 2, 4, 8, ... seconds;
    once `max_attempts` is reached the failure is logged and the object is
    returned without a record. No exception reaches the caller.

    Args:
        instance: The object whose accession is fetched.
        online_database: The Entrez database to fetch from, e.g. "protein".
        _attempt: The number of the current attempt. Defaults to 1.
        max_attempts: The number of attempts before giving up. Defaults to
            the configured max_retrieval_attempts (3 when unset).
        _entrez_email: The email address sent to Entrez. Defaults to the
            configured entrez_email.

    Returns:
        The same object, with its GenBank record set when the fetch
        succeeded and unchanged otherwise. See [incomplete_dict_cleaner] in
        utils.py to remove objects left without a record.
    """
    Entrez.email = _entrez_email  # type: ignore[assignment]

    kwargs = {
        "db": online_database,
        "id": str(instance.accession),
        "rettype": "gb",
        "retmode": "text",
    }
    if instance.HSP:
        kwargs |= {
            "seq_start": max(1, instance.HSP.sbjct_start),
            "seq_stop": instance.HSP.sbjct_end,
        }
    try:
        with Entrez.efetch(**kwargs) as handle:  # type: ignore[no-untyped-call]
            genbank_record = handle.read()
            instance.set_genbank(genbank_record)
        return instance

    except Exception as e:
        if _attempt < max_attempts:
            time.sleep(2**_attempt)
            # Retries are routine with NCBI; only the final failure matters.
            logger.debug(
                f"{instance.accession}: GenBank fetch failed ({e!s}); retry {_attempt + 1}"
            )
            return gb_fetcher(instance, online_database, _attempt + 1, max_attempts)
        logger.error(
            f"{instance.accession}: GenBank record not fetched after {max_attempts} "
            f"attempts ({e!s}). The probe will lack its sequence."
        )
        return instance


def gb_executor(
    object_dict: dict[str, RetroSeeker],
    online_database: str,
    max_attempts: int = defaults.MAX_RETRIEVAL_ATTEMPTS,
) -> dict[str, RetroSeeker] | None:
    """Fetch the GenBank record of every object in a dictionary, one at a time.

    Args:
        object_dict: The objects to fetch records for.
        online_database: The Entrez database to fetch from, e.g. "protein".
        max_attempts: The number of attempts per object before giving up.
            Defaults to the configured max_retrieval_attempts.

    Returns:
        The same objects keyed by ``{accession}-{identifier}``, each with its
        GenBank record when the fetch succeeded, or None (after a CRITICAL log
        line) when `object_dict` is empty.
    """
    full_retrieved_results = {}

    with tqdm(
        total=len(object_dict), desc="Fetching GenBank sequences", disable=None
    ) as object_bar:
        for value in object_dict.values():
            # gb_fetcher always returns the instance (updated on success, or
            # unchanged after exhausting retries) - never None.
            result = gb_fetcher(
                instance=value,
                online_database=online_database,
                max_attempts=max_attempts,
            )
            key_identifier = f"{value.accession}-{value.identifier}"
            full_retrieved_results[key_identifier] = result
            logger.debug(
                f"Added {key_identifier} to GenBank Dictionary\n{result.display_info()}"
            )
            object_bar.update(1)

    if not full_retrieved_results:
        logger.critical("No fetched GenBank results. Exiting.")
        return None

    return full_retrieved_results
