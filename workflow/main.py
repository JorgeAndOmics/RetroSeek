import subprocess
import os

import defaults

import logging
from colored_logging import colored_logging
def main():
    logging.debug('Starting workflow.')
    logging.debug('Starting BLAST Database Generation')
    subprocess.run(['python', 'db_generator.py'])
    # logging.debug('Starting LTR Retrieval')
    # subprocess.run(['python', 'ltr_retriever.py'])
    # logging.debug('Starting Probe Extraction')
    # subprocess.run(['python', 'probe_extractor.py'])
    # logging.debug('Starting Species BLAST Sequence Retrieval')
    # subprocess.run(['python', 'species_seq_retriever.py'])
    logging.debug('Starting Virus BLAST Sequence Retrieval')
    subprocess.run(['python', 'virus_seq_retriever.py'])
    logging.debug('Starting Reciprocal BLAST Database Generation')
    subprocess.run(['python', 'rec_db_generator.py'])
    logging.debug('Starting Reciprocal BLAST')
    subprocess.run(['python', 'rec_blaster.py'])
    logging.debug('Starting Domain Retrieval')
    subprocess.run(['python', 'domain_retriever.py'])

if __name__ == '__main__':
    colored_logging('main.txt')
    main()