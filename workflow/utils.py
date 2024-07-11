from ratelimit import limits, sleep_and_retry
import pickle
import random
import string
import os

import defaults

import logging


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
        pickle.dump(data, f)


def unpickler(input_directory_path, input_file_name: str) -> dict:
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
            :raises FileNotFoundError: If the file does not exist.
            :raises pickle.UnpicklingError: If the file is not a valid pickle file or is corrupted.
            :raises Exception: If an error occurs while unpickling the file.

    """
    try:
        file_path = os.path.join(input_directory_path, input_file_name)
        with open(file_path, 'rb') as f:
            return pickle.load(f)
    except FileNotFoundError as e:
        raise FileNotFoundError(f"The file {file_path} does not exist.") from e
    except pickle.UnpicklingError:
        raise pickle.UnpicklingError(f"The file {file_path} is not a valid pickle file or is corrupted.")
    except Exception as e:
        raise Exception(f"An error occurred while unpickling the file {file_path}: {e}")


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


@sleep_and_retry
@limits(calls=defaults.MAX_EXECUTION_ATTEMPTS_PER_SECOND, period=defaults.MIN_EXECUTION_INTERVAL)
def execution_limiter(func, *args, **kwargs):
    """
    Limits the number of executions of a function per second.

        Parameters
        ----------
            :param func: The function to limit.

    """
    return func(*args, **kwargs)


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
