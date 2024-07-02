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

    probe_dict = {
        f'{str(row["Accession"])}': Object(
            family=str(row['Family']),
            virus=str(row['Name']),
            abbreviation=str(row['Abbreviation']),
            probe=str(row['Probe']),
            accession=str(row['Accession']),
            identifier=random_string_generator(6)
        )
        for index, row in probe_table.iterrows()
    }

    logging.info(f'Parsed {len(probe_dict)} probes')
    return probe_dict


def probe_extractor(probe_dict: dict, online_database: str) -> dict:
    """
    Orchestrates the probe extraction process, retrieves GenBank information and pickles
    the dictionary containing the probe objects.

        Parameters
        ----------
            :param probe_dict: The dictionary containing the probe objects.
            :param online_database: The online database to be used for sequence retrieval.

        Returns
        -------
            :returns: None

        Raises
        ------
            :raises Exception: If the retrieval process fails.

    """
    for value in probe_dict.values():
        try:
            gb_fetcher(value, online_database)
            logging.info(f'Successfully retrieved {value.accession}-{value.identifier} GenBank information')
        except Exception as e:
            logging.warning(f'Failed to retrieve {value.accession}-{value.identifier} GenBank information: {e}')

    return probe_dict


if __name__ == '__main__':
    colored_logging(log_file_name='probe_extractor.txt')

    probe_dict: dict = table_parser(input_csv_file=os.path.join('..', 'data', 'tables', 'Probes.csv'))

    probe_extraction: dict = probe_extractor(probe_dict=probe_dict,
                                             online_database='protein')

    incomplete_dict_cleaner(object_dict=probe_extraction)

    pickler(data=probe_extraction,
            output_directory_path=defaults.PICKLE_DIR,
            output_file_name='probe_dict.pkl')

    logging.info('Probe extraction and retrieval completed')
