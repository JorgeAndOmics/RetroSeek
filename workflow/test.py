import os
import io
from io import StringIO
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import db_utils
import pprint

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


dict_test = print(defaults.REVERSE_FULL_DICT)