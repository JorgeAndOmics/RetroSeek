"""
This script sets up colored logging for the application.

It configures the logging module to output logs to both a file and the console,
with colored output for easier readability. The log file name is specified as a parameter.

Modules:
    - coloredlogs: Provides colored logging for the console.
    - logging: Standard Python logging module.
    - os: Provides a way of using operating system dependent functionality.
    - defaults: Custom module for default configurations.

Functions:
    - colored_logging(log_file_name: str) -> None: Configures logging with colored output.
"""

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
