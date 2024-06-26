import os
import io
from io import StringIO
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import defaults
from utils import pickler, unpickler, colored_logging
from object_class import Object
from seq_utils import blaster, blaster_parser, seq_fetcher

from Bio.Blast import NCBIXML
from Bio import Entrez

import logging, coloredlogs

test = unpickler(os.path.join('..', 'data', 'pickles'), 'probe_dict.pkl')
# print(test['AJG42161.1'].get_gff())

handle = test['AJG42161.1'].get_fasta()

print(handle)