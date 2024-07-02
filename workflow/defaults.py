import os

# Genomes
SPECIES: list = ['Rhinolophus_ferrumequinum']
VIRUS: list = ['NCBI_virus']

# BLAST
EXPANSION_SIZE: int = 0
E_VALUE: float = 0.1
ACCESSION_ID_REGEX = r'[A-Z]{2,}_?\d*\.\d{1,2}'
PROBE_MIN_LENGTH: dict = {
    'GAG': 400,
    'POL': 200,
    'ENV': 400,
    'VIG': 400,
    'N_protein': 400,
    'P_protein': 400,
    'G_protein': 400,
    'L_protein': 400,
    'X_protein': 400,
    'M_protein': 400
}

# Logging
LEVEL_STYLES = {
    'debug': {'color': 'white'},     # Standard debug level
    'info': {'color': 'cyan', 'bold': 'yes'},     # Standard info level
    'warning': {'color': 'yellow'},  # Standard warning level
    'error': {'color': 'red', 'bold': 'yes'},      # Standard error level
    'critical': {'color': 'black', 'bold': 'yes', 'background': 'red'},  # Standard critical level
}

FIELD_STYLES = {
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
ROOT_DB = os.path.join('/', 'mnt', 'v', 'databases', 'local')
SPECIES_DB = os.path.join(ROOT_DB, 'species')
VIRUS_DB = os.path.join(ROOT_DB, 'virus')

# Directories
LOG_DIR = os.path.join('..', 'logs')
PICKLE_DIR = os.path.join('..', 'data', 'pickles')
TMP_DIR = os.path.join('..', 'data', 'tmp')

# Execution and requests
MAX_RETRIEVAL_ATTEMPTS = 9
MAX_EXECUTION_ATTEMPTS_PER_SECOND = 10
MIN_EXECUTION_INTERVAL = 1  # seconds
ENTREZ_EMAIL: str = 'jgonzlez@tcd.ie'
NCBI_API_TOKEN = 'faa9e17bb461e82963f079c167ec5c7aac08'

# Displays
DISPLAY_REQUESTS_WARNING = False
