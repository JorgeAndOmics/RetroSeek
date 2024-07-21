import os
import io
from io import StringIO
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import db_utils

import defaults
from utils import pickler, unpickler
from object_class import Object
from colored_logging import colored_logging

from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq

import logging, coloredlogs


probe_test = unpickler(input_directory_path=defaults.PICKLE_DIR,
                       input_file_name='probe_dict.pkl')

db_utils.objdict2fasta(object_dict=probe_test,
                       output_directory_path=defaults.PICKLE_DIR,
                       output_filename='probe_test.fasta')