import os
import sys
import shutil

SPECIES: list = ['Rhinolophus_ferrumequinum']
PROBES: list = ['gag',
                'pol',
                'env',
                'vif',
                'n_protein',
                'p_protein',
                'g_protein',
                'l_protein',
                'x_protein',
                'm_protein']
ENTREZ_EMAIL: str = 'jgonzlez@tcd.ie'
EXPANSION_SIZE: int = 5000
E_VALUE: float = 0.1
EXPANSION_SWITCH: str = 'N'
ACCESSION_ID_REGEX = '[A-Z]*\d*[._]\d*'
PROBE_MIN_LENGTH: dict = {
    'gag': 0,
    'pol': 0,
    'env': 0,
    'vif': 0,
    'n_protein': 0,
    'p_protein': 0,
    'g_protein': 0,
    'l_protein': 0,
    'x_protein': 0,
    'm_protein': 0
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
FASTA_DIR = os.path.join('data', 'fastas')
LOG_DIR = os.path.join('..', 'logs')
PICKLE_DIR = os.path.join('..', 'pickles')
TMP_DIR = os.path.join('..', 'data', 'tmp')
