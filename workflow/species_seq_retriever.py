import os
import io
from io import StringIO
import subprocess
from collections import defaultdict

import defaults
from utils import pickler, unpickler, directory_content_eraser
from object_class import Object
from seq_utils import *
from colored_logging import colored_logging

from Bio.Blast import NCBIXML
from Bio import Entrez

import logging, coloredlogs



if __name__ == '__main__':
    colored_logging('species_seq_retriever.txt')

    probe_dict = unpickler(os.path.join('..', 'data', 'pickles'), 'probe_dict.pkl')

    blast2gb2pickle(object_dict=probe_dict,
                    command='tblastn',
                    online_database='nucleotide',
                    input_database_path=defaults.SPECIES_DB,
                    species=defaults.SPECIES,
                    output_pickle_directory_path=defaults.PICKLE_DIR,
                    output_pickle_file_name='tblastn_results.pkl')

    directory_content_eraser(defaults.TMP_DIR)
