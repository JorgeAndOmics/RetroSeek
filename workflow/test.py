import os
import io
from io import StringIO
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

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


test = unpickler(os.path.join('..', 'data', 'pickles'), 'tblastn_results.pkl')

print(len(test.items()))