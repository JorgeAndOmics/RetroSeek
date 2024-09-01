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

dict = unpickler(defaults.PICKLE_DIR, 'tblastn_results.pkl')

full_dict = {}
species_dict = {}

for key, value in dict.items():
    species_dict[value.species] = {key: value}

species_list = list(set(species_dict.keys()))

for key, value in species_dict():  # key is species name, value is dictionary of objects with that species
    species_dict = seq_utils.blast_retriever(object_dict=value,  # value is the dictionary of objects with that species
                                             command='blastn',
                                             genome=[value.species],
                                             online_database='nucleotide',
                                             input_database_path=defaults.LTR_DB,  # LTR_DB is the directory where the LTRHarvest output is stored
                                             multi_threading=True)


full_dict |= species_dict

pickler(full_dict, defaults.PICKLE_DIR, 'ltrharvest_blast_against.pkl')