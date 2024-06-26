import pickle
import os

import defaults
from object_class import Object


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




