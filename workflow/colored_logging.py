import coloredlogs
import logging
import os

import defaults


def colored_logging(log_file_name: str) -> None:
    """
    Sets up logging and configures coloredlogs with the custom fields and level styles

        Parameters
        ----------
            :param log_file_name: The name of the file to save the log in.
    """
    # Configure coloredlogs with the custom field and level styles
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(message)s', handlers=[
        logging.FileHandler(os.path.join(defaults.LOG_DIR, log_file_name), mode='w'),
        logging.StreamHandler()
    ]
                        )

    coloredlogs.install(
        level='DEBUG',
        fmt='%(asctime)s - %(message)s',
        level_styles=defaults.LEVEL_STYLES,
        field_styles=defaults.FIELD_STYLES
    )
