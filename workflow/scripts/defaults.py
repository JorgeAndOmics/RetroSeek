"""
Defaults Configuration Script
=============================

This script loads configuration settings from a YAML file and sets up various
constants and directory paths used throughout the project.

Modules:
    - `yaml`: For parsing the YAML configuration file.
    - `os`: For handling file and directory paths.

Configuration:
    - The configuration file is expected to be located at `../data/config/config.yaml`.
    - Various constants and directory paths are initialized based on the configuration file.

Usage:
    This script is intended to be imported as a module and not run directly.
"""

import yaml
import os

CONFIG_FILE = os.path.abspath(os.path.join('..', 'data', 'config', 'config.yaml'))
with open(CONFIG_FILE, 'r') as f:
    config = yaml.safe_load(f)

# BLAST
E_VALUE = config['blast']['e_value']
ACCESSION_ID_REGEX = r'[A-Z]{2,}_?[0-9]+\.[0-9]{1,2}'
PROBE_MIN_LENGTH = config['parameters']['probe_min_length']

# Logging
LEVEL_STYLES = config['logging']['level_styles']
FIELD_STYLES = config['logging']['field_styles']

# Directories
PATH_DICT = {
    'ROOT': os.path.abspath(config['root'].get('db_root_folder', os.path.join('..')))
}

# === Root Directories ===
PATH_DICT['DATA_DIR'] = os.path.abspath(os.path.join(config['root'].get('data_root_folder', os.path.join(PATH_DICT['ROOT'], 'data'))))
PATH_DICT['RESULTS_DIR'] = os.path.abspath(os.path.join(config['root'].get('results_root_folder', os.path.join(PATH_DICT['ROOT'], 'results'))))
PATH_DICT['LOG_DIR'] = os.path.abspath(os.path.join(config['root'].get('logs_root_folder', os.path.join(PATH_DICT['ROOT'], 'logs'))))
PATH_DICT['WORKFLOW_DIR'] = os.path.abspath(os.path.join(config['root'].get('workflow_root_folder', os.path.join(PATH_DICT['ROOT'], 'workflow'))))

# === Workflow Directories ===
PATH_DICT['SCRIPTS_DIR'] = os.path.abspath(os.path.join(PATH_DICT['WORKFLOW_DIR'], 'scripts'))

# === Database Directories ===
PATH_DICT['ROOT_DB'] = os.path.abspath(os.path.join(PATH_DICT['ROOT'], 'local'))
PATH_DICT['SPECIES_DB'] = os.path.abspath(os.path.join(PATH_DICT['ROOT_DB'], 'blast_dbs', 'species'))
PATH_DICT['ACCESSORY_DB'] = os.path.abspath(os.path.join(PATH_DICT['ROOT'], 'accessory'))

# === Data Subdirectories ===
PATH_DICT['CONFIG_DIR'] = os.path.abspath(os.path.join(PATH_DICT['DATA_DIR'], 'config'))
PATH_DICT['SPECIES_DIR'] = os.path.abspath(os.path.join(PATH_DICT['DATA_DIR'], 'species'))
PATH_DICT['TABLE_INPUT_DIR'] = os.path.abspath(os.path.join(PATH_DICT['DATA_DIR'], 'tables'))
PATH_DICT['PICKLE_DIR'] = os.path.abspath(os.path.join(PATH_DICT['DATA_DIR'], 'pickles'))
PATH_DICT['TMP_DIR'] = os.path.abspath(os.path.join(PATH_DICT['DATA_DIR'], 'tmp'))
PATH_DICT['TBLASTN_PICKLE_DIR'] = os.path.abspath(os.path.join(PATH_DICT['PICKLE_DIR'], 'tblastn'))

# === Results - Tables ===
PATH_DICT['TABLE_OUTPUT_DIR'] = os.path.abspath(os.path.join(PATH_DICT['RESULTS_DIR'], 'tables'))
PATH_DICT['TABLE_HOTSPOT_DIR'] = os.path.abspath(os.path.join(PATH_DICT['TABLE_OUTPUT_DIR'], 'hotspots'))
PATH_DICT['TABLE_OVERLAP_MATRIX_DIR'] = os.path.abspath(os.path.join(PATH_DICT['TABLE_OUTPUT_DIR'], 'overlap_matrix'))
PATH_DICT['SEGMENTED_SPECIES_DIR'] = os.path.abspath(os.path.join(PATH_DICT['TABLE_OUTPUT_DIR'], 'segmented_species'))
PATH_DICT['PLOT_DATAFRAMES_DIR'] = os.path.abspath(os.path.join(PATH_DICT['TABLE_OUTPUT_DIR'], 'plot_dataframes'))

# === Results - Plots ===
PATH_DICT['PLOT_DIR'] = os.path.abspath(os.path.join(PATH_DICT['RESULTS_DIR'], 'plots'))
PATH_DICT['CIRCLE_PLOT_DIR'] = os.path.abspath(os.path.join(PATH_DICT['PLOT_DIR'], 'circle_plots'))
PATH_DICT['HOTSPOT_PDF_DIR'] = os.path.abspath(os.path.join(PATH_DICT['PLOT_DIR'], 'hotspot_pdfs'))

# === Results - Tracks ===
PATH_DICT['TRACK_DIR'] = os.path.abspath(os.path.join(PATH_DICT['RESULTS_DIR'], 'tracks'))
PATH_DICT['TRACK_ORIGINAL_DIR'] = os.path.abspath(os.path.join(PATH_DICT['TRACK_DIR'], 'original'))
PATH_DICT['TRACK_CANDIDATES_DIR'] = os.path.abspath(os.path.join(PATH_DICT['TRACK_DIR'], 'candidates'))
PATH_DICT['TRACK_VALID_DIR'] = os.path.abspath(os.path.join(PATH_DICT['TRACK_DIR'], 'valid'))
PATH_DICT['TRACK_HOTSPOTS_DIR'] = os.path.abspath(os.path.join(PATH_DICT['TRACK_DIR'], 'hotspots'))

# === Results - LTR ===
PATH_DICT['LTRHARVEST_DIR'] = os.path.abspath(os.path.join(PATH_DICT['RESULTS_DIR'], 'ltrharvest'))
PATH_DICT['LTRDIGEST_DIR'] = os.path.abspath(os.path.join(PATH_DICT['RESULTS_DIR'], 'ltrdigest'))

# === Logs & Workflow ===
PATH_DICT['DOWNLOAD_LOG'] = os.path.abspath(os.path.join(PATH_DICT['LOG_DIR'], 'download_log.log'))

# === Accessory Tools ===
PATH_DICT['HMM_PROFILE_DIR'] = os.path.abspath(os.path.join(PATH_DICT['ACCESSORY_DB'], 'hmm_profiles'))


# Directory generation
for value in PATH_DICT.values():
    os.makedirs(value, exist_ok=True)

# Execution and requests
NUM_CORES = config['execution'].get('num_cores', 1)
USE_SPECIES_DICT = config.get('execution', False).get('use_species_dict', False)
RETRIVAL_TIME_LAG = config['execution'].get('retrieval_time_lag', 0.3)
MAX_RETRIEVAL_ATTEMPTS = config['execution'].get('max_retrieval_attempts', 3)
MAX_THREADPOOL_WORKERS = config['execution'].get('max_threadpool_workers', 1)
ENTREZ_EMAIL = config['execution'].get('entrez_email', '')

# Display
DISPLAY_SNAKEMAKE_INFO: bool = config['display']['display_snakemake_info']
DISPLAY_REQUESTS_WARNING: bool = config['display']['display_requests_warning']
DISPLAY_OPERATION_INFO: bool = config['display']['display_operation_info']

# INPUT
PROBE_CSV = os.path.abspath(config['input'].get('probe_csv'))

# Genomes
SPECIES_DICT: dict = config.get('species', {})

if not USE_SPECIES_DICT:
    SPECIES: list = [(f.split('.fa')[0]).strip() for f in os.listdir(SPECIES_DB) if f.endswith('.fa')]
else:
    SPECIES: list = SPECIES_DICT.keys()
