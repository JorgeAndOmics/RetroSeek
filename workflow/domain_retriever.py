import utils
import defaults
import online_seq_utils

from colored_logging import colored_logging
import logging


if __name__ == '__main__':
    colored_logging(log_file_name='domain_retriever.txt')

    virus_dict: dict = utils.unpickler(input_directory_path=defaults.PICKLE_DIR,
                                       input_file_name='virus_blastn_results.pkl')

    clean_domain_results: dict = online_seq_utils.online_blast_threadpool_executor(object_dict=virus_dict,
                                                                                   command='rpsblast',
                                                                                   online_database='cdd')

    utils.pickler(data=clean_domain_results,
                  output_directory_path=defaults.PICKLE_DIR,
                  output_file_name='domain_results.pkl')

    utils.directory_content_eraser(directory_path=defaults.TMP_DIR)

    logging.info(f'Successfully performed Virus BLASTn and retrieval on {len(clean_domain_results)} sequences.')
