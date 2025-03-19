from ratelimit import limits, sleep_and_retry
import cloudpickle
import pickle
import random
import string
import dill
import os

from functools import wraps
from tqdm import tqdm

import defaults

import logging

import Bio

def pickler(data, output_directory_path, output_file_name: str) -> None:
    """
    Pickles the data to the specified path.

        Parameters
        ----------
            :param data: The data to be pickled.
            :param output_directory_path: The directory where the file will be saved.
            :param output_file_name: The name of the file to save the pickled data in.

    """
    if not os.path.exists(os.path.join(output_directory_path)):
        os.makedirs(os.path.join(output_directory_path), exist_ok=True)
    with open(os.path.join(output_directory_path, output_file_name), 'wb') as f:
        dill.dump(data, f)


def unpickler(input_directory_path, input_file_name: str):
    """
    Unpickles the data from the specified path.

        Parameters
        ----------
            :param input_directory_path: The directory where the file is saved.
            :param input_file_name: The name of the file to load the pickled data from.

        Returns
        -------
            :returns: The unpickled data.
            
        Raises
        ------
            :raises Exception: If the unpickling process fails.

    """
    try:
        file_path = os.path.join(input_directory_path, input_file_name)
        with open(file_path, 'rb') as f:
            return dill.load(f)
    except Exception as e:
        logging.error(f'Failed to unpickle {file_path}')
        raise Exception(f'Failed to unpickle {file_path}') from e


def directory_generator(parent_directory_path,
                        new_directory_name):
    """
    Generates a directory if it does not exist.

        Parameters
        ----------
            :param parent_directory_path: The directory in which to generate the new subdirectory.
            :param new_directory_name: The name of the new subdirectory.

        Returns
        -------
            :returns: The path to the generated directory

    """
    new_path = os.path.join(parent_directory_path, new_directory_name)
    if not os.path.exists(new_path):
        os.makedirs(new_path, exist_ok=True)
    return new_path


def directory_file_retriever(input_directory_path) -> list:
    """
    Retrieve all files in a directory (not directories).

    Args:
        :param input_directory_path: Path to the directory containing the files.

    Returns:
        :returns: List of files in the directory.
    """
    return [f for f in os.listdir(input_directory_path) if os.path.isfile(os.path.join(input_directory_path, f))]


def directory_content_eraser(directory_path) -> None:
    """
    Erases the content of the specified directory.

    Parameters
    ----------
        :param directory_path: The directory to erase the content of.

    """
    logging.debug('Cleaning up temporal files...')

    for file in os.listdir(directory_path):
        file_path = os.path.join(directory_path, file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
        except Exception as e:
            logging.warning(f'Failed to delete {file_path}: {e}')

def get_last_directory(input_directory_path):
    """
    Returns the name of the last directory in the input path.
        
        Parameters
        ----------
            :param input_directory_path: The path to the directory.
            
    """
    # Ensure the path ends with a separator to handle the case where the path is a directory
    if not input_directory_path.endswith(os.path.sep):
        path = os.path.join(input_directory_path, '')

    # Get the directory name of the path (remove the trailing separator)
    parent_path = os.path.dirname(input_directory_path.rstrip(os.path.sep))

    return os.path.basename(parent_path)

def incomplete_dict_cleaner(object_dict: dict) -> dict:
    """
    Removes incomplete objects from an object dictionary.

    Parameters
    ----------
        :param object_dict: The dictionary to remove incomplete records from.

    Returns
    -------
        :returns: The input dictionary without incomplete records.

    """
    logging.debug('Cleaning up incomplete objects...')

    return {key: value for key, value in object_dict.items() if value.is_complete()}

def random_string_generator(length: int) -> str:
    """
    Generates a random string of the specified length.

        Parameters
        ----------
        :param length: The length of the string to generate.

        Returns
        -------
        :returns: The generated random string.

    """
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=length))



















































































































































































































































