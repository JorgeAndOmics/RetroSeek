from ratelimit import limits, sleep_and_retry
import pickle
import random
import string
import os
import re

import defaults
from object_class import Object

import logging, coloredlogs

def pickler(data, output_directory_path, output_file_name)-> None:
    """
    Pickles the data to the specified path.

    Args:
        data (Any): The data to be pickled.
        output_directory_path (str): The directory where the file will be saved.
        output_file_name (str): The name of the file to save the pickled data in.

    Raises:
        OSError: If the directory cannot be created.
        IOError: If the file cannot be written.
    """
    if not os.path.exists(os.path.join(output_directory_path)):
        os.makedirs(os.path.join(output_directory_path), exist_ok=True)
    with open(os.path.join(output_directory_path, output_file_name), 'wb') as f:
        pickle.dump(data, f)

def unpickler(input_directory_path, input_file_name):
    """
    Unpickles the data from the specified path.

    Args:
        input_directory_path (str): The directory where the file is saved.
        input_file_name (str): The name of the file to load the pickled data from.

    Returns:
        Any: The unpickled data.

    Raises:
        FileNotFoundError: If the file cannot be found.
        IOError: If the file cannot be read.
    """
    with open(os.path.join(input_directory_path, input_file_name), 'rb') as f:
        return pickle.load(f)

def directory_content_eraser(directory_path)-> None:
    """
    Erases the content of the specified directory.

    Args:
        directory_path (str): The directory to erase the content of.

    Raises:
        OSError: If the directory cannot be created.
        IOError: If the file cannot be written.
    """
    logging.debug('Cleaning up temporal files...')

    for file in os.listdir(directory_path):
        file_path = os.path.join(directory_path, file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except Exception as e:
            logging.warning(f'Failed to delete {file_path}: {e}')


def incomplete_dict_cleaner(object_dict: dict)-> dict:
    """
    Removes incomplete objects from an object dictionary.

    Args:
        object_dict (dict): The dictionary to remove incomplete records from.

    Returns:
        dict: The input dictionary without incomplete records.
    """
    logging.debug('Cleaning up incomplete objects...')

    return {key: value for key, value in object_dict.items() if value.is_complete()}

@sleep_and_retry
@limits(calls=defaults.MAX_EXECUTION_ATTEMPTS_PER_SECOND, period=defaults.MIN_EXECUTION_INTERVAL)
def execution_limiter(func,*args, **kwargs)-> None:
    """
    Limits the number of executions of a function per second.

        Args:
            func (function): The function to limit.

        Returns:
            Any: The result of the function.
    """
    return func(*args, **kwargs)

def random_string_generator(length: int)-> str:
    """
    Generates a random string of the specified length.

    Args:
        length (int): The length of the string to generate.

    Returns:
        str: The generated string.
    """
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=length))

def accession_finder(query:str, pattern:str)-> str:
    """
    Finds the accession ID in a query using a regex pattern.

    Args:
        query (str): The query to search the accession ID in.
        pattern (str): The regex pattern to use for the search.

    Returns:
        str: The accession ID found in the query.

    CAUTION!: The function returns the first occurrence matching the pattern.
    """
    regex_pattern = re.compile(pattern)
    return str(regex_pattern.search(query).group())



