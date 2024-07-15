import pandas as pd
import numpy as np
import os
import logging

import defaults
import utils
from colored_logging import colored_logging


def df_generator(column_names) -> pd.DataFrame:
    """
    Generates an empty dataframe with the specified columns.

            Returns
            -------
                :returns: A DataFrame with the specified columns.

    """
    return pd.DataFrame(columns=column_names)


def df_initializer(attributes: list, input_dataframe: pd.DataFrame, input_dictionary: dict) -> pd.DataFrame:
    """
    Initializes a pre-generated DataFrame with the input dictionary data.

        Parameters
        ----------
            :param attributes: The list of object attributes to initialize the DataFrame with. It's the same as the
            columns in the DataFrame, but with the first letter capitalized.
            :param input_dataframe: The DataFrame to initialize.
            :param input_dictionary: The object pair dictionary to initialize the DataFrame with.

        Returns
        -------
            :returns: A DataFrame with the data from the object pair dictionary.

    """
    for key, objct in input_dictionary.items():
        for attribute in attributes:
            input_dataframe.loc[key, attribute] = getattr(objct, attribute.lower())
            logging.info(f'Added {key} -> {attribute} to DataFrame')

    return input_dataframe


if __name__ == '__main__':

    files = 'probe_dict'

    colored_logging(log_file_name='obj2dict.txt')

    objct_dict: dict = utils.unpickler(input_directory_path=defaults.PICKLE_DIR,
                                       input_file_name=f'{files}.pkl')

    df: pd.DataFrame = df_generator(column_names=defaults.CSV_ATTRIBUTES)

    df: pd.DataFrame = df_initializer(attributes=defaults.CSV_ATTRIBUTES,
                                      input_dataframe=df,
                                      input_dictionary=objct_dict)

    logging.info(f'Generated DataFrame with {len(df)} rows.')
    print(df)

    df.to_csv(os.path.join(defaults.TABLE_OUTPUT_DIR, f'{files}.csv'))
