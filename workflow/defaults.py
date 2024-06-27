import os
import sys
import shutil

SPECIES: list = ['Rhinolophus_ferrumequinum']
ENTREZ_EMAIL: str = 'jgonzlez@tcd.ie'
EXPANSION_SIZE: int = 0
E_VALUE: float = 0.1
ACCESSION_ID_REGEX = '[A-Z]{2,}_?\d*\.\d{1,2}'
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

# Define a custom level style to color INFO level logs green
LEVEL_STYLES = {
    'debug': {'color': 'white'},     # Standard debug level
    'info': {'color': 'cyan', 'bold': 'yes'},     # Standard info level
    'warning': {'color': 'yellow'}, # Standard warning level
    'error': {'color': 'red', 'bold': 'yes'},      # Standard error level
    'critical': {'color': 'black', 'bold': 'yes', 'background': 'red'},  # Standard critical level
}

# Define custom styles to color fields green and levels with their respective colors
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

SPECIES_DB = '/mnt/v/databases/local'
VIRUS_DB = '/mnt/v/databases/ncbi_virus_db/ncbi_virus'
VIRUS_FASTA = '/mnt/v/databases/refseq_virus_db/viral.1.1.genomic.fna'
LOG_DIR = os.path.join('..', 'logs')
PICKLE_DIR = os.path.join('..', 'data', 'pickles')
TMP_DIR = os.path.join('..', 'data', 'tmp')
MAX_RETRIEVAL_ATTEMPTS = 7