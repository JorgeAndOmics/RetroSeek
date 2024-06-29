from seq_utils import *
from utils import *

from colored_logging import colored_logging
import logging

if __name__ == '__main__':
    colored_logging(log_file_name='virus_seq_retriever.txt')

    tblastn_results = unpickler(input_directory_path=os.path.join('..', 'data', 'pickles'),
                           input_file_name='tblastn_results.pkl')

    clean_blastn2gb_results = blast_retriever(object_dict=tblastn_results,
                                              command='blastn',
                                              genome=defaults.VIRUS,
                                              online_database='nucleotide',
                                              input_database_path=defaults.VIRUS_DB,
                                              multi_threading=True)

    pickler(data=clean_blastn2gb_results,
            output_directory_path=defaults.PICKLE_DIR,
            output_file_name='virus_blastn_results.pkl')

    directory_content_eraser(directory_path=defaults.TMP_DIR)

    logging.info(f'Successfully performed Virus BLASTn and retrieval on {len(clean_tblastn2gb_results)} sequences.')