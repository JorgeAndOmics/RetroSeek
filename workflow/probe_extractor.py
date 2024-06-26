import os
import pandas as pd

from utils import pickler, unpickler, colored_logging
from object import Object
from seq_utils import seq_fetcher

import logging, coloredlogs

def table_parser(input_csv_file):
    '''
    Parse the CSV files containing the probe genes and their respective IDs
    '''
    probe_dict = {}
    # Read the CSV file
    probe_table = pd.read_csv(input_csv_file)
    for index, row in probe_table.iterrows():
        probe_dict[str(row['Accession'])] = Object(
                                  family=str(row['Family']),
                                  virus=str(row['Name']),
                                  abbreviation=str(row['Abbreviation']),
                                  probe=str(row['Probe']),
                                  accession=str(row['Accession'])
                                              )

    logging.info(f'Parsed {len(probe_dict)} probes')
    return probe_dict

def probe_extractor(input_csv_file, output_pickle_directory_path, output_pickle_file_name, online_database):
    '''
    Orchestrates the probe extraction process
    '''
    probe_dict = table_parser(input_csv_file)
    for key, value in probe_dict.items():
        seq_fetcher(value, online_database)
    pickler(probe_dict, output_pickle_directory_path, output_pickle_file_name)

if __name__ == '__main__':
    colored_logging('probe_extractor.txt')

    probe_extractor(os.path.join('..', 'data', 'tables', 'Probes.csv'),
                    os.path.join('..', 'data', 'pickles'),
                    'probe_dict.pkl',
                    'protein')

    logging.info(f'Probe extraction completed')





