"""
Object-to-Table Converter
=========================

This script deserializes one or more pickled Object dictionaries, extracts their attributes,
and compiles them into species-specific DataFrames, which are then concatenated and saved
as both CSV and Parquet files.

Modules:
    - `utils`: Utility functions.
    - `defaults`: project-wide constants such as thread count and directory paths.
    - `colored_logging`: for structured logging.
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
    python obj2table.py --files genomeA.pkl genomeB.pkl --output_file output/probes
"""

# =============================================================================
# Imports and Logging Setup
# =============================================================================
import pandas as pd
import os
import logging
import argparse
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

import defaults
import utils
from colored_logging import colored_logging

# =============================================================================
# 1. Attribute Extraction Function
# =============================================================================
def extract_attributes_from_object(obj) -> dict:
    """
    Extracts all structured attributes from a probe Object for tabular export.

        Parameters
        ----------
            :param obj: A single Object instance representing a probe. Must have `.genbank`, `.alignment`, and `.HSP` as optional attributes.

        Returns
        -------
            :returns: data: A dictionary mapping field names to scalar values, ready for DataFrame construction. It includes:
                - Basic metadata from the Object itself
                - GenBank fields (if `.genbank` is not None)
                - Alignment metadata (if `.alignment` is not None)
                - BLAST HSP features (if `.HSP` is not None)

        Raises
        ------
            :raises AttributeError: If the Object structure does not expose required attributes.
            :raises Exception: For any unexpected failure during data extraction.
    """
    data: dict = {
        # Basic Object attributes
        'label': obj.label,
        'virus': obj.virus,
        'abbreviation': obj.abbreviation,
        'species': obj.species,
        'probe': obj.probe,
        'accession': obj.accession,
        'identifier': obj.identifier,
        'strand': obj.strand,
    }
    if defaults.USE_SPECIES_DICT:
        data |= {
            'species_name': defaults.SPECIES_DICT.get(obj.species, obj.species),
        }
    else:
        data |= {
            'species_name': obj.species,
        }

    # GenBank info
    gb = obj.genbank
    data |= {
        'genbank_id': gb.id if gb else None,
        'genbank_name': gb.name if gb else None,
        'genbank_description': gb.description if gb else None,
        'genbank_dbxrefs': gb.dbxrefs if gb else None,
        'genbank_annotations': str(gb.annotations) if gb else None,
        'genbank_seq': str(gb.seq) if gb else None,
    }

    # Alignment info
    aln = obj.alignment
    data |= {
        'alignment_title': aln.title if aln else None,
        'alignment_length': aln.length if aln else None,
        'alignment_accession': aln.accession if aln else None,
        'alignment_hit_id': aln.hit_id if aln else None,
        'alignment_hit_def': aln.hit_def if aln else None,
    }

    # HSP info
    hsp = obj.HSP
    data |= {
        'hsp_bits': hsp.bits if hsp else None,
        'hsp_score': hsp.score if hsp else None,
        'hsp_evalue': hsp.expect if hsp else None,
        'hsp_query': hsp.query if hsp else None,
        'hsp_sbjct': hsp.sbjct if hsp else None,
        'hsp_query_start': hsp.query_start if hsp else None,
        'hsp_query_end': hsp.query_end if hsp else None,
        'hsp_sbjct_start': hsp.sbjct_start if hsp else None,
        'hsp_sbjct_end': hsp.sbjct_end if hsp else None,
        'hsp_identity': hsp.identities if hsp else None,
        'hsp_align_length': hsp.align_length if hsp else None,
        'hsp_gaps': hsp.gaps if hsp else None,
        'hsp_positives': hsp.positives if hsp else None,
        'hsp_strand': hsp.strand if hsp else None,
        'hsp_frame': hsp.frame if hsp else None,
    }

    return data


# =============================================================================
# 2. Main Execution Block
# =============================================================================
if __name__ == '__main__':
    # Initialize logging
    colored_logging(log_file_name='obj2dict.txt')

    # -------------------------------------------------------------------------
    # 2.1 Argument Parsing
    # -------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='Converts serialized probe Objects into a tabular DataFrame.')
    parser.add_argument(
        '--files',
        type=str,
        nargs='+',
        required=True,
        help='One or more pickle files to load from PICKLE_DIR.'
    )
    parser.add_argument(
        '--output_file_name',
        type=str,
        required=True,
        help='Base name for output files (CSV and Parquet will be generated).'
    )
    args = parser.parse_args()

    # -------------------------------------------------------------------------
    # 2.2 Load Object Dictionaries from Pickles
    # -------------------------------------------------------------------------
    all_objects: list = []
    for file in args.files:
        objct_dict = utils.unpickler(
            input_directory_path=defaults.PATH_DICT['PICKLE_DIR'],
            input_file_name=file
        )
        logging.info(f'{os.path.basename(file).split(".")[0]}: {len(objct_dict)} objects retrieved')
        all_objects.extend(objct_dict.values())

    logging.info(f'Total objects loaded: {len(all_objects)}')

    # -------------------------------------------------------------------------
    # 2.3 Group Objects by Species
    # -------------------------------------------------------------------------
    species_objects: dict = defaultdict(list)
    for obj in all_objects:
        species_objects[obj.species].append(obj)

    logging.info(f'Objects grouped into {len(species_objects)} species')

    # -------------------------------------------------------------------------
    # 2.4 Process Each Species Sequentially and Build DataFrames
    # -------------------------------------------------------------------------
    species_dataframes: dict = {}

    for species, objects in species_objects.items():
        logging.info(f'Processing species: {species} ({len(objects)} objects)')

        results: list[dict] = []
        with tqdm(total=len(objects), desc=f'Processing {species}') as pbar:
            with ThreadPoolExecutor(max_workers=defaults.MAX_THREADPOOL_WORKERS) as executor:
                futures = [executor.submit(extract_attributes_from_object, obj) for obj in objects]
                for future in as_completed(futures):
                    results.append(future.result())
                    pbar.update()

        species_df = pd.DataFrame(results)
        species_dataframes[species] = species_df
        logging.info(f'Species {species}: DataFrame created with {len(species_df)} rows')

    # -------------------------------------------------------------------------
    # 2.5 Concatenate All Species DataFrames
    # -------------------------------------------------------------------------
    logging.info('Concatenating all species DataFrames...')
    df = pd.concat(species_dataframes.values(), ignore_index=True)
    logging.info(f'Final DataFrame: {len(df)} total rows from {len(species_dataframes)} species')

    # -------------------------------------------------------------------------
    # 2.6 Save DataFrame to CSV and Parquet
    # -------------------------------------------------------------------------
    output_csv_path = f'{args.output_file_name}.csv'
    output_parquet_path = f'{args.output_file_name}.parquet'

    df.to_csv(output_csv_path, index=False)
    logging.info(f'CSV saved to: {output_csv_path}')

    df.to_parquet(output_parquet_path, index=False)
    logging.info(f'Parquet saved to: {output_parquet_path}')