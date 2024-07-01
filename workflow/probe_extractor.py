from colored_logging import colored_logging
from object_class import Object
from seq_utils import *
from utils import pickler
import pandas as pd
import defaults
import logging
import os


def table_parser(input_csv_file):
    """
    Parse the CSV files containing the probe genes and their respective IDs

        Parameters
        ----------
            :param input_csv_file: The path to the CSV file containing probe info.

        Returns
        -------
            :returns: A dictionary containing object pairs.
    """
    # Read the CSV file
    probe_table = pd.read_csv(input_csv_file)

    random_string = random_string_generator(6)

    probe_dict = {
            f'{str(row["Accession"])}-{random_string}': Object(
            family=str(row['Family']),
            virus=str(row['Name']),
            abbreviation=str(row['Abbreviation']),
            probe=str(row['Probe']),
            accession=str(row['Accession']),
            identifier=random_string
        )
        for index, row in probe_table.iterrows()
    }
    logging.info(f'Parsed {len(probe_dict)} probes')
    return probe_dict


def probe_extractor(input_csv_file, output_pickle_directory_path, output_pickle_file_name, online_database):
    """
    Orchestrates the probe extraction process, retrieves GenBank information and pickles
    the dictionary containing the probe objects.

        Parameters
        ----------
            :param input_csv_file: The path to the CSV file containing probe info.
            :param output_pickle_directory_path: The path to the directory where the pickled file will be saved.
            :param output_pickle_file_name: The name of the pickled file to be saved.
            :param online_database: The online database to be used for sequence retrieval.

        Returns
        -------
            :returns: None
    """
    probe_dict = table_parser(input_csv_file)
    for key, value in probe_dict.items():
        gb_fetcher(value, online_database)
    pickler(probe_dict, output_pickle_directory_path, output_pickle_file_name)

if __name__ == '__main__':
    colored_logging(log_file_name='probe_extractor.txt')

    probe_extractor(input_csv_file=os.path.join('..', 'data', 'tables', 'Probes.csv'),
                    output_pickle_directory_path=defaults.PICKLE_DIR,
                    output_pickle_file_name='probe_dict.pkl',
                    online_database='protein')

    logging.info('Probe extraction completed')





