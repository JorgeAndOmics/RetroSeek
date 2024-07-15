import os

# Genomes
SPECIES: list = ['Rhinolophus_ferrumequinum']
VIRUS: list = ['NCBI_virus']

# BLAST
EXPANSION_SIZE: int = 0
E_VALUE: float = 0.09
ACCESSION_ID_REGEX: str = r'[A-Z]{2,}_?\d*\.\d{1,2}'
PROBE_MIN_LENGTH: dict = {
    'GAG': 200,
    'POL': 200,
    'ENV': 200,
    'VIF': 200,
    'N_protein': 200,
    'P_protein': 200,
    'G_protein': 200,
    'L_protein': 200,
    'X_protein': 200,
    'M_protein': 200
}

# Logging
LEVEL_STYLES: dict = {
    'debug': {'color': 'white'},     # Standard debug level
    'info': {'color': 'cyan', 'bold': 'yes'},     # Standard info level
    'warning': {'color': 'yellow'},  # Standard warning level
    'error': {'color': 'red', 'bold': 'yes'},      # Standard error level
    'critical': {'color': 'black', 'bold': 'yes', 'background': 'red'},  # Standard critical level
}

FIELD_STYLES: dict = {
    'asctime': {'color': 'green'},
    'hostname': {'color': 'green'},
    'levelname': {'color': 'green'},
    'name': {'color': 'green'},
    'programname': {'color': 'green'},
    'username': {'color': 'green'},
    'process': {'color': 'green'},
    'thread': {'color': 'green'}
}

# Databases
ROOT_DB = os.path.join('/mnt/', 'v', 'databases', 'local')
SPECIES_DB = os.path.join(ROOT_DB, 'species')
VIRUS_DB = os.path.join(ROOT_DB, 'virus')

# Directories
TABLE_INPUT_DIR = os.path.join('..', 'data', 'tables', 'input')
TABLE_OUTPUT_DIR = os.path.join('..', 'data', 'tables', 'output')
LOG_DIR = os.path.join('..', 'logs')
PICKLE_DIR = os.path.join('..', 'data', 'pickles')
TMP_DIR = os.path.join('..', 'data', 'tmp')

# Execution and requests
MAX_RETRIEVAL_ATTEMPTS: int = 9
MAX_EXECUTION_ATTEMPTS_PER_SECOND: int = 10
MIN_EXECUTION_INTERVAL: int = 1  # seconds
ENTREZ_EMAIL: str = 'jgonzlez@tcd.ie'
NCBI_API_TOKEN: str = 'faa9e17bb461e82963f079c167ec5c7aac08'

# Displays
DISPLAY_REQUESTS_WARNING: bool = False

# CSV
CSV_DELIMITER: str = ','
CSV_ATTRIBUTES: list[str] = ['Family',
                             'Virus',
                             'Abbreviation',
                             'Species',
                             'Probe',
                             'Accession',
                             'Identifier',
                             'Strand']
