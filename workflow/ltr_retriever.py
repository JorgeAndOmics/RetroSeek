"""
This script generates the Suffixerator index, and uses it to generate the LTRHarvest output for the input directory.
"""
import os

import defaults
import utils
import db_utils

from colored_logging import colored_logging


if __name__ == '__main__':
    # Set up logging
    colored_logging('ltr_retriever.txt')

    # TODO: Finish implementation: Follow the trace of db_generator, from files to directories.

