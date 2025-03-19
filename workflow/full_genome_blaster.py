import seq_utils
import utils
import defaults

import argparse

from colored_logging import colored_logging
import logging

if __name__ == '__main__':
    colored_logging(log_file_name='full_genome_blaster.txt')

    # Add argument parser
    parser = argparse.ArgumentParser(description='Performs tBLASTn on a genome.')
    parser.add_argument('--genome', type=str, help='The name of the genome to perform tBLASTn on.')
    parser.add_argument('--num_threads', type=str, default=defaults.MAX_THREADPOOL_WORKERS, help='The number of threads to use for tBLASTn.')
    args = parser.parse_args()
    genome = args.genome
    num_threads = args.num_threads

    probe_dict: dict = utils.unpickler(input_directory_path=defaults.PICKLE_DIR,
                                       input_file_name='probe_dict.pkl')

    tblastn_results: dict = seq_utils.blast_retriever(object_dict=probe_dict,
                                                               command='tblastn',
                                                               genome=args.genome,
                                                               input_database_path=defaults.SPECIES_DB,
                                                               num_threads=num_threads,
                                                               display_full_info=False)

    utils.pickler(data=tblastn_results,
                  output_directory_path=defaults.TBLASTN_PICKLE_DIR,
                  output_file_name=f'{genome}.pkl')

    # logging.info(f'Successfully performed tBLASTn and retrieval on {genome}: {len(tblastn_results)} sequences retrieved.')