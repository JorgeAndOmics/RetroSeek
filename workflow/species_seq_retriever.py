from seq_utils import *
from utils import *

from colored_logging import colored_logging
import logging

if __name__ == '__main__':
    colored_logging(log_file_name='species_seq_retriever.txt')

    probe_dict = unpickler(input_directory_path=os.path.join('..', 'data', 'pickles'),
                           input_file_name='probe_dict.pkl')

    clean_tblastn2gb_results = blast_retriever(object_dict=probe_dict,
                                               command='tblastn',
                                               genome=defaults.SPECIES,
                                               online_database='nucleotide',
                                               input_database_path=defaults.SPECIES_DB,
                                               multi_threading=True)

    pickler(data=clean_tblastn2gb_results,
            output_directory_path=defaults.PICKLE_DIR,
            output_file_name='tblastn_results.pkl')

    directory_content_eraser(directory_path=defaults.TMP_DIR)

    logging.info(f'Successfully performed tBLASTn and retrieval on {len(clean_tblastn2gb_results)} sequences.')