import os
import io
from io import StringIO
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import db_utils
import pprint
import cloudpickle
import pickle
import seq_utils

import defaults
from utils import pickler, unpickler
from object_class import Object
from colored_logging import colored_logging

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Blast import Record

import logging, coloredlogs

import sys
import utils
from Bio import Entrez

from Bio.Blast import NCBIWWW, NCBIXML

test_dict = utils.unpickler(input_directory_path=defaults.PICKLE_DIR, input_file_name='full_genome_blast.pkl')

print(len(test_dict.keys()))
