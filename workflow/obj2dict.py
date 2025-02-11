import pandas as pd
import os
import logging

import defaults
import utils
from colored_logging import colored_logging

from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import argparse


def extract_attributes_from_object(obj):
    """
    Extracts attributes from the Object instance.
    If tabular_data is available, use it to override fields that can be derived from tabular results.
    """
    data = {}
    # Extract attributes from the Object itself
    data['family'] = obj.family
    data['virus'] = obj.virus
    data['abbreviation'] = obj.abbreviation
    data['species'] = obj.species
    data['probe'] = obj.probe
    data['accession'] = obj.accession
    data['identifier'] = obj.identifier
    data['strand'] = obj.strand

    # Extract attributes from genbank (SeqRecord)
    if obj.genbank:
        gb = obj.genbank
        data['genbank_id'] = gb.id
        data['genbank_name'] = gb.name
        data['genbank_description'] = gb.description
        data['genbank_dbxrefs'] = gb.dbxrefs
        data['genbank_annotations'] = str(gb.annotations)
        data['genbank_seq'] = str(gb.seq)
    else:
        data['genbank_id'] = None
        data['genbank_name'] = None
        data['genbank_description'] = None
        data['genbank_dbxrefs'] = None
        data['genbank_annotations'] = None
        data['genbank_seq'] = None

    # Extract attributes from alignment (from XML)
    if obj.alignment:
        aln = obj.alignment
        data['alignment_title'] = aln.title
        data['alignment_length'] = aln.length
        data['alignment_accession'] = aln.accession
        data['alignment_hit_id'] = aln.hit_id
        data['alignment_hit_def'] = aln.hit_def
    else:
        data['alignment_title'] = None
        data['alignment_length'] = None
        data['alignment_accession'] = None
        data['alignment_hit_id'] = None
        data['alignment_hit_def'] = None

    # Extract attributes from HSP (from XML)
    if obj.HSP:
        hsp = obj.HSP
        data['hsp_bits'] = hsp.bits
        data['hsp_score'] = hsp.score
        data['hsp_evalue'] = hsp.expect
        data['hsp_query'] = hsp.query
        data['hsp_sbjct'] = hsp.sbjct
        data['hsp_query_start'] = hsp.query_start
        data['hsp_query_end'] = hsp.query_end
        data['hsp_sbjct_start'] = hsp.sbjct_start
        data['hsp_sbjct_end'] = hsp.sbjct_end
        data['hsp_identity'] = hsp.identities
        data['hsp_align_length'] = hsp.align_length
        data['hsp_gaps'] = hsp.gaps
        data['hsp_positives'] = hsp.positives
        data['hsp_strand'] = hsp.strand
        data['hsp_frame'] = hsp.frame
    else:
        data['hsp_bits'] = None
        data['hsp_score'] = None
        data['hsp_evalue'] = None
        data['hsp_query'] = None
        data['hsp_sbjct'] = None
        data['hsp_query_start'] = None
        data['hsp_query_end'] = None
        data['hsp_sbjct_start'] = None
        data['hsp_sbjct_end'] = None
        data['hsp_identity'] = None
        data['hsp_align_length'] = None
        data['hsp_positives'] = None
        data['hsp_gaps'] = None
        data['hsp_strand'] = None
        data['hsp_frame'] = None

    return data


if __name__ == '__main__':
    colored_logging(log_file_name='obj2dict.txt')

    parser = argparse.ArgumentParser(description='Converts objects to DataFrame.')
    parser.add_argument('--file', type=str, default=None, help='Name of the pickle file in PKL_DIR to load.')
    args = parser.parse_args()

    if args.file:
        files = args.file
        filename = files.split(".")[0]
    else:
        logging.error('No file name provided. Exiting script.')
        exit(1)

    objct_dict: dict = utils.unpickler(
        input_directory_path=defaults.PICKLE_DIR,
        input_file_name=f'{files}'
    )

    logging.info(f'Loaded {len(objct_dict)} objects from pickle file.')

    objects = list(objct_dict.values())

    # Number of parallel workers
    max_workers_num = defaults.MAX_THREADPOOL_WORKERS

    results = []
    with tqdm(total=len(objects), desc='Processing objects') as pbar:
        with ProcessPoolExecutor(max_workers=max_workers_num) as executor:
            futures = [executor.submit(extract_attributes_from_object, obj) for obj in objects]
            for future in as_completed(futures):
                results.append(future.result())
                pbar.update()

    df = pd.DataFrame(results)

    output_csv_path = os.path.join(defaults.TABLE_OUTPUT_DIR, f'{filename}.csv')
    output_parquet_path = os.path.join(defaults.TABLE_OUTPUT_DIR, f'{filename}.parquet')

    df.to_csv(output_csv_path, index=False)
    logging.info(f'Saved DataFrame as CSV to {output_csv_path}')

    df.to_parquet(output_parquet_path, index=False)
    logging.info(f'Saved DataFrame as Parquet to {output_parquet_path}')
