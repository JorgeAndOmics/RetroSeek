import pickle
import os

import defaults
from object_class import Object


def pickler(data, output_directory_path, output_file_name):
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

def directory_content_eraser(directory_path):
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


def incomplete_dict_cleaner(object_dict: dict):
    """
    Removes incomplete objects from an object dictionary.

    Args:
        object_dict (dict): The dictionary to clean.

    Returns:
        dict: The cleaned dictionary.
    """
    logging.debug('Cleaning up incomplete objects...')

    return {key: value for key, value in object_dict.items() if value.is_complete()}



