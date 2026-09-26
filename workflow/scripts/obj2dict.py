r"""Object-to-table converter: pickled hit objects to one CSV and Parquet table.

This script deserializes one or more pickled Object dictionaries, extracts their attributes,
and compiles them into species-specific DataFrames, which are then concatenated and saved
as both CSV and Parquet files.

Modules:
    - `utils`: Utility functions.
    - `defaults`: project-wide constants such as thread count and directory paths.
    - `log`: the console line contract and job log (ADR-021).
    - `argparse`, `pandas`, `os`, `tqdm`, `concurrent.futures`: for CLI, tabular processing, I/O, and parallelism.

Main Workflow:
    1. Parse command-line arguments for input pickle files and output base name.
    2. Unpickle and load Object dictionaries.
    3. Group objects by species.
    4. Extract attributes from each object (including nested GenBank, Alignment, and HSP data).
    5. Store results in species-specific DataFrames within a dictionary.
    6. Concatenate all species DataFrames into a unified DataFrame.
    7. Save the final DataFrame to both CSV and Parquet formats.

Requirements:
    - Input `.pkl` files must reside in `defaults.PICKLE_DIR`.
    - The `Object` class must be readable; no changes in environment or class definition between pickling and unpickling.

Usage:
    python obj2dict.py --files genomeA.pkl genomeB.pkl \
        --csv_path results/tables/<name>/<name>.csv \
        --parquet_path data/tables/<name>/<name>.parquet
"""

# =============================================================================
# Imports and Logging Setup
# =============================================================================
import argparse
import logging
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd
from tqdm import tqdm

import defaults
import utils
from log import OK, job_logging, run_main

logger = logging.getLogger(__name__)


# =============================================================================
# 1. Attribute Extraction Function
# =============================================================================
# Column name -> the attribute it is read from, for the alignment and HSP records
# attached to a probe Object. A missing record gives None in every one of its
# columns. Plain attribute names rather than getter functions: this runs once
# per BLAST hit, hundreds of thousands of times per genome.
_ALIGNMENT_COLUMNS = {
    "alignment_title": "title",
    "alignment_length": "length",
    "alignment_accession": "accession",
    "alignment_hit_id": "hit_id",
    "alignment_hit_def": "hit_def",
}
_HSP_COLUMNS = {
    "hsp_bits": "bits",
    "hsp_score": "score",
    "hsp_evalue": "expect",
    "hsp_query": "query",
    "hsp_sbjct": "sbjct",
    "hsp_query_start": "query_start",
    "hsp_query_end": "query_end",
    "hsp_sbjct_start": "sbjct_start",
    "hsp_sbjct_end": "sbjct_end",
    "hsp_identity": "identities",
    "hsp_align_length": "align_length",
    "hsp_gaps": "gaps",
    "hsp_positives": "positives",
    "hsp_strand": "strand",
    "hsp_frame": "frame",
}


def _read_columns(record: Any, columns: dict[str, str]) -> dict[str, Any]:
    """One optional record's columns; all None when the record is absent."""
    if not record:
        return dict.fromkeys(columns)
    return {name: getattr(record, attr) for name, attr in columns.items()}


def _genbank_columns(gb: Any) -> dict[str, Any]:
    """The GenBank record's columns; annotations and sequence as text."""
    if not gb:
        return dict.fromkeys(
            [
                "genbank_id",
                "genbank_name",
                "genbank_description",
                "genbank_dbxrefs",
                "genbank_annotations",
                "genbank_seq",
            ]
        )
    return {
        "genbank_id": gb.id,
        "genbank_name": gb.name,
        "genbank_description": gb.description,
        "genbank_dbxrefs": gb.dbxrefs,
        "genbank_annotations": str(gb.annotations),
        "genbank_seq": str(gb.seq),
    }


def _species_name(species: str) -> str:
    """The configured display name of a genome, or its stem without a species map."""
    if defaults.USE_SPECIES_DICT:
        return defaults.SPECIES_DICT.get(species, species)
    return species


def extract_attributes_from_object(obj: Any) -> dict[str, Any]:
    """Flatten one probe object into a table row for export.

    Args:
        obj: One RetroSeeker object. Its `genbank`, `alignment` and `HSP` may
            each be unset.

    Returns:
        A dictionary mapping column names to values, ready for DataFrame
        construction. Every row has the same columns: the object's own
        metadata (with the accession cut to its first word) and species
        display name, then the GenBank, alignment and HSP columns, which are
        None when that record is unset.

    Raises:
        AttributeError: If the object lacks an attribute a column is read from.
    """
    # Defensive accession cleanup: older pickles (pre-seq_utils-fix) stored
    # the full BLAST hit_def ("CM138268.1 Molossus molossus chromosome 3,
    # ...") as accession. Split on whitespace and keep the first token -
    # idempotent for already-clean accessions. Keeps the parquet seqid
    # column compatible with LTRdigest / GRanges seqnames downstream.
    raw_accession = obj.accession or ""
    clean_accession = raw_accession.split()[0] if raw_accession else raw_accession
    return {
        "label": obj.label,
        "virus": obj.virus,
        "abbreviation": obj.abbreviation,
        "species": obj.species,
        "probe": obj.probe,
        "accession": clean_accession,
        "identifier": obj.identifier,
        "strand": obj.strand,
        "species_name": _species_name(obj.species),
        **_genbank_columns(obj.genbank),
        **_read_columns(obj.alignment, _ALIGNMENT_COLUMNS),
        **_read_columns(obj.HSP, _HSP_COLUMNS),
    }


# =============================================================================
# 2. Main Execution Block
# =============================================================================
def main(args: argparse.Namespace) -> None:
    """Load the pickled objects and write them as one CSV and one Parquet table."""
    # -------------------------------------------------------------------------
    # 2.2 Load Object Dictionaries from Pickles
    # -------------------------------------------------------------------------
    all_objects: list[Any] = []
    for file in args.files:
        objct_dict = utils.unpickler(
            input_directory_path=defaults.PATH_DICT["PICKLE_DIR"], input_file_name=file
        )
        logger.info(f"{Path(file).stem}: {len(objct_dict)} objects retrieved")
        all_objects.extend(objct_dict.values())

    logger.info(f"Total objects loaded: {len(all_objects)}")

    # -------------------------------------------------------------------------
    # 2.3 Group Objects by Species
    # -------------------------------------------------------------------------
    species_objects: dict[str, list[Any]] = defaultdict(list)
    for obj in all_objects:
        species_objects[obj.species].append(obj)

    logger.info(f"Objects grouped into {len(species_objects)} species")

    # -------------------------------------------------------------------------
    # 2.4 Process Each Species Sequentially and Build DataFrames
    # -------------------------------------------------------------------------
    species_dataframes: dict[str, Any] = {}

    for species, objects in species_objects.items():
        logger.info(f"Processing species: {species} ({len(objects)} objects)")

        results: list[dict[str, Any]] = []
        with (
            tqdm(
                total=len(objects), desc=f"Processing {species}", disable=None
            ) as pbar,
            ThreadPoolExecutor(max_workers=defaults.MAX_THREADPOOL_WORKERS) as executor,
        ):
            futures = [
                executor.submit(extract_attributes_from_object, obj) for obj in objects
            ]
            for future in as_completed(futures):
                results.append(future.result())
                pbar.update()

        species_df = pd.DataFrame(results)
        species_dataframes[species] = species_df
        logger.info(f"Species {species}: DataFrame created with {len(species_df)} rows")

    # -------------------------------------------------------------------------
    # 2.5 Concatenate All Species DataFrames
    # -------------------------------------------------------------------------
    logger.info("Concatenating all species DataFrames...")
    df = pd.concat(species_dataframes.values(), ignore_index=True)
    logger.info(
        f"Final DataFrame: {len(df)} total rows from {len(species_dataframes)} species"
    )

    # -------------------------------------------------------------------------
    # 2.6 Save DataFrame to CSV and Parquet
    # -------------------------------------------------------------------------
    df.to_csv(args.csv_path, index=False)
    logger.info(f"CSV saved to: {args.csv_path}")

    df.to_parquet(args.parquet_path, index=False)
    logger.info(f"Parquet saved to: {args.parquet_path}")
    logger.log(
        OK,
        "%s objects from %s species tabulated",
        f"{len(df):,}",
        len(species_dataframes),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Converts serialized probe Objects into a tabular DataFrame."
    )
    parser.add_argument(
        "--files",
        type=str,
        nargs="+",
        required=True,
        help="One or more pickle files to load from PICKLE_DIR.",
    )
    parser.add_argument(
        "--csv_path",
        type=str,
        required=True,
        help="Full path for the user-facing CSV output.",
    )
    parser.add_argument(
        "--parquet_path",
        type=str,
        required=True,
        help="Full path for the pipeline-internal Parquet output.",
    )
    parser.add_argument("--log", type=Path, help="job log (the Snakemake log: path)")
    args = parser.parse_args()
    # Its two rules (probe_extractor, and the blast_pkl2parquet checkpoint) must
    # stay byte-identical, so neither passes --log: rerunning either would wake
    # heavy work. The log follows the usual layout under the script's name.
    job_log = args.log or Path(defaults.PATH_DICT["LOG_DIR"]) / "obj2dict" / "all.log"
    job_logging(job_log, "obj2dict")
    run_main(lambda: main(args))
