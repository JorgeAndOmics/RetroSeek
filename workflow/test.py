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


def extract_strand_from_HSP(HSP_object) -> str:
    """
    Extracts the strand information from the HSP object

        Returns:
            None

        Raises:
            None
    """
    try:
        if HSP_object.frame:
            if HSP_object.frame[-1] > 0 and isinstance(HSP_object.frame[-1], int):
                return '+'
            if HSP_object.frame[-1] < 0 and isinstance(HSP_object.frame[-1], int):
                return '-'
        else:
            if HSP_object.sbjct_start < HSP_object.sbjct_end:
                return '+'
            if HSP_object.sbjct_start > HSP_object.sbjct_end:
                return '-'
    except Exception as e:
        logging.warning(f'Could not extract strand: {e}')
        return None

test = unpickler(os.path.join('..', 'data', 'pickles'), 'tblastn_results.pkl')
# print(test['AJG42161.1'].get_gff())

test_instance = (test['CM014239.1-CEPMGO'])

test_hsp = test_instance.HSP

test_frame = test_hsp.frame

test_frame_method = extract_strand_from_HSP(test_hsp)

print(test_hsp)
