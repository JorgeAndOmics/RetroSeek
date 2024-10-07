import pandas as pd
import os
import logging

import defaults
import utils
from colored_logging import colored_logging

from multiprocessing import Pool, Queue, cpu_count
from tqdm import tqdm


def extract_attributes_from_object(obj):
    """
    Extracts attributes from the Object instance and returns a dictionary.
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

    # Extract attributes from alignment
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

    # Extract attributes from HSP
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

    files = 'ltr_fasta_blast'

    colored_logging(log_file_name='obj2dict.txt')

    logging.info('Starting obj2dict script.')

    objct_dict: dict = utils.unpickler(input_directory_path=defaults.PICKLE_DIR,
                                       input_file_name=f'{files}.pkl')

    logging.info(f'Loaded {len(objct_dict)} objects from pickle file.')

    # Prepare to process objects in parallel
    objects = list(objct_dict.values())

    # Use multiprocessing Pool
    num_workers = 1
    logging.info(f'Using {num_workers} parallel workers for processing.')

    def init_pool(tqdm_instance):
        global pbar
        pbar = tqdm_instance

    # TODO: TRANSFORM TQDM CONTEXT INTO WRAPPER FUNCTION/DECORATOR
    with tqdm(total=len(objects), desc='Processing objects') as pbar:
        with Pool(processes=num_workers, initializer=init_pool, initargs=(pbar,)) as pool:
            results = []
            for res in pool.imap_unordered(extract_attributes_from_object, objects):
                results.append(res)
                pbar.update()
            pool.close()
            pool.join()

    logging.info('Finished processing all objects.')

    # Create DataFrame
    df = pd.DataFrame(results)

    logging.info(f'Generated DataFrame with {len(df)} rows.')

    output_csv_path = os.path.join(defaults.TABLE_OUTPUT_DIR, f'{files}.csv')
    output_parquet_path = os.path.join(defaults.TABLE_OUTPUT_DIR, f'{files}.parquet')

    df.to_csv(output_csv_path, index=False)
    logging.info(f'Saved DataFrame as CSV to {output_csv_path}')

    df.to_parquet(output_parquet_path, index=False)
    logging.info(f'Saved DataFrame as Parquet to {output_parquet_path}')

    logging.info('Script completed successfully.')
