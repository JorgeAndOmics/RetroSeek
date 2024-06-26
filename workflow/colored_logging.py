import os

import defaults

import logging, coloredlogs



def colored_logging(file_name: str):
    """
    Sets up logging and configures coloredlogs with the custom fields and level styles

        Args:
            file_name: The name of the file to save the log in.

        Returns:
            None
    """
    # Configure coloredlogs with the custom field and level styles
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(message)s', handlers=[
        logging.FileHandler(os.path.join(defaults.LOG_DIR, file_name), mode='w'),
        logging.StreamHandler()
    ]
                        )

    coloredlogs.install(
        level='DEBUG',
        fmt='%(asctime)s - %(message)s',
        level_styles=defaults.LEVEL_STYLES,
        field_styles=defaults.FIELD_STYLES
    )