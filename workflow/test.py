import os
import io
from io import StringIO
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import defaults
from utils import pickler, unpickler
from object_class import Object
from seq_utils import blaster, blaster_parser, gb_fetcher
from colored_logging import colored_logging

from Bio.Blast import NCBIXML
from Bio import Entrez

import logging, coloredlogs



test = unpickler(os.path.join('..', 'data', 'pickles'), 'tblastn_results.pkl')
# print(test['AJG42161.1'].get_gff())

for value in test.values():
    print(value.strand)