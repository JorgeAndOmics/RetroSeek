import pandas as pd
import os
from io import StringIO
import time
import subprocess
import tempfile
from dataclasses import dataclass, field
from typing import Optional, Any
import defaults
import pickle
import re

import random
import string

from Bio import SeqIO, Entrez
from Bio.Blast import NCBIXML, NCBIWWW

import logging, coloredlogs

def colored_logging(file_name: str):
    """
    Sets up logging and configures coloredlogs with the custom fields and level styles
    
        Args:
            file_name: The name of the file to save the log in.
        
        Returns:
            None
    """
    # Configure coloredlogs with the custom field and level styles
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(message)s', handlers=[
        logging.FileHandler(os.path.join(defaults.LOG_DIR, file_name), mode='w'),
        logging.StreamHandler()
    ]
                        )

    coloredlogs.install(
        level='DEBUG',
        fmt='%(asctime)s - %(message)s',
        level_styles=defaults.LEVEL_STYLES,
        field_styles=defaults.FIELD_STYLES
    )

def pickler(data, directory_path, file_name):
    """
    Pickles the data to the specified path.

    Args:
        data (Any): The data to be pickled.
        directory_path (str): The directory where the file will be saved.
        file_name (str): The name of the file to save the pickled data in.

    Raises:
        OSError: If the directory cannot be created.
        IOError: If the file cannot be written.
    """
    if not os.path.exists(os.path.join(directory_path)):
        os.makedirs(os.path.join(directory_path), exist_ok=True)
    with open(os.path.join(directory_path, file_name), 'wb') as f:
        pickle.dump(data, f)

def unpickler(directory_path, file_name):
    """
    Unpickles the data from the specified path.

    Args:
        directory_path (str): The directory where the file is saved.
        file_name (str): The name of the file to load the pickled data from.

    Returns:
        Any: The unpickled data.

    Raises:
        FileNotFoundError: If the file cannot be found.
        IOError: If the file cannot be read.
    """
    with open(os.path.join(directory_path, file_name), 'rb') as f:
        return pickle.load(f)




