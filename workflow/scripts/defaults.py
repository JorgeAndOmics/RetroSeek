"""
Defaults Configuration Script
=============================

This script loads configuration settings from a YAML file and sets up various
constants and directory paths used throughout the project.

Modules:
    - `yaml`: For parsing the YAML configuration file.
    - `pathlib.Path`: For handling file and directory paths.

Configuration:
    - The configuration file is expected to be located at `../data/config/config.yaml`.
    - Various constants and directory paths are initialized based on the configuration file.

Usage:
    This script is intended to be imported as a module and not run directly.
"""
from pathlib import Path
import yaml

CONFIG_FILE = Path(__file__).parents[2] / 'data' / 'config' / 'config.yaml'

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
    'ROOT': Path(config['root'].get('db_root_folder', '..')).resolve()
}

# === Root Directories ===
PATH_DICT['DATA_DIR'] = Path(config['root'].get('data_root_folder', PATH_DICT['ROOT'] / 'data')).resolve()
PATH_DICT['RESULTS_DIR'] = Path(config['root'].get('results_root_folder', PATH_DICT['ROOT'] / 'results')).resolve()
PATH_DICT['LOG_DIR'] = Path(config['root'].get('logs_root_folder', PATH_DICT['ROOT'] / 'logs')).resolve()
PATH_DICT['WORKFLOW_DIR'] = Path(__file__).parent

# === Workflow Directories ===
PATH_DICT['SCRIPTS_DIR'] = (PATH_DICT['WORKFLOW_DIR'] / 'scripts').resolve()

# === Database Directories ===
PATH_DICT['ROOT_DB'] = (PATH_DICT['ROOT']).resolve()
PATH_DICT['SPECIES_DB'] = (PATH_DICT['ROOT_DB']).resolve()
PATH_DICT['ACCESSORY_DB'] = (PATH_DICT['ROOT'] / 'accessory').resolve()

# === Data Subdirectories ===
PATH_DICT['CONFIG_DIR'] = (PATH_DICT['DATA_DIR'] / 'config').resolve()
PATH_DICT['SPECIES_DIR'] = (PATH_DICT['DATA_DIR'] / 'species').resolve()
PATH_DICT['TABLE_INPUT_DIR'] = (PATH_DICT['DATA_DIR'] / 'tables').resolve()
PATH_DICT['PICKLE_DIR'] = (PATH_DICT['DATA_DIR'] / 'pickles').resolve()
PATH_DICT['TMP_DIR'] = (PATH_DICT['DATA_DIR'] / 'tmp').resolve()
PATH_DICT['TBLASTN_PICKLE_DIR'] = (PATH_DICT['PICKLE_DIR'] / 'tblastn').resolve()

# === Results - Tables ===
PATH_DICT['TABLE_OUTPUT_DIR'] = (PATH_DICT['RESULTS_DIR'] / 'tables').resolve()
PATH_DICT['TABLE_HOTSPOT_DIR'] = (PATH_DICT['TABLE_OUTPUT_DIR'] / 'hotspots').resolve()
PATH_DICT['TABLE_OVERLAP_MATRIX_DIR'] = (PATH_DICT['TABLE_OUTPUT_DIR'] / 'overlap_matrix').resolve()
PATH_DICT['SEGMENTED_SPECIES_DIR'] = (PATH_DICT['TABLE_OUTPUT_DIR'] / 'segmented_species').resolve()
PATH_DICT['PLOT_DATAFRAMES_DIR'] = (PATH_DICT['TABLE_OUTPUT_DIR'] / 'plot_dataframes').resolve()
PATH_DICT['TABLE_PAIR_DIR'] = (PATH_DICT['TABLE_OUTPUT_DIR'] / 'probe_pairs').resolve()

# === Results - Plots ===
PATH_DICT['PLOT_DIR'] = (PATH_DICT['RESULTS_DIR'] / 'plots').resolve()
PATH_DICT['CIRCLE_PLOT_DIR'] = (PATH_DICT['PLOT_DIR'] / 'circle_plots').resolve()
PATH_DICT['HOTSPOT_PDF_DIR'] = (PATH_DICT['PLOT_DIR'] / 'hotspot_pdfs').resolve()

# === Results - Tracks ===
PATH_DICT['TRACK_DIR'] = (PATH_DICT['RESULTS_DIR'] / 'tracks').resolve()
PATH_DICT['TRACK_ORIGINAL_DIR'] = (PATH_DICT['TRACK_DIR'] / 'original').resolve()
PATH_DICT['TRACK_CANDIDATES_DIR'] = (PATH_DICT['TRACK_DIR'] / 'candidates').resolve()
PATH_DICT['TRACK_VALID_DIR'] = (PATH_DICT['TRACK_DIR'] / 'valid').resolve()
PATH_DICT['TRACK_HOTSPOTS_DIR'] = (PATH_DICT['TRACK_DIR'] / 'hotspots').resolve()

# === Results - LTR ===
PATH_DICT['LTRHARVEST_DIR'] = (PATH_DICT['TRACK_DIR'] / 'ltrharvest').resolve()
PATH_DICT['LTRDIGEST_DIR'] = (PATH_DICT['TRACK_DIR'] / 'ltrdigest').resolve()
PATH_DICT['SOLO_LTR_DIR'] = (PATH_DICT['TRACK_DIR'] / 'solo_ltr').resolve()
PATH_DICT['FLANKING_LTR_DIR'] = (PATH_DICT['TRACK_DIR'] / 'flanking_ltr').resolve()

# === Logs & Workflow ===
PATH_DICT['DOWNLOAD_LOG'] = (PATH_DICT['LOG_DIR'] / 'download_log.log').resolve()

# === Accessory Tools ===
PATH_DICT['HMM_PROFILE_DIR'] = (PATH_DICT['ACCESSORY_DB'] / 'hmm_profiles').resolve()


# Directory generation
for value in PATH_DICT.values():
    value.mkdir(parents=True, exist_ok=True)

# Execution and requests
NUM_CORES = config['execution'].get('num_cores', 1)
USE_SPECIES_DICT = config.get('execution', False).get('use_species_dict', False)
RETRIVAL_TIME_LAG = config['execution'].get('retrieval_time_lag', 0.3)
MAX_RETRIEVAL_ATTEMPTS = config['execution'].get('max_retrieval_attempts', 3)
MAX_THREADPOOL_WORKERS = config['execution'].get('max_threadpool_workers', 1)
ENTREZ_EMAIL = config['execution'].get('entrez_email', '')

# Display
DISPLAY_SNAKEMAKE_INFO: bool = config['display'].get('display_snakemake_info', False)
DISPLAY_REQUESTS_WARNING: bool = config['display'].get('display_requests_warning', False)
DISPLAY_OPERATION_INFO: bool = config['display'].get('display_operation_info', False)

# INPUT
PROBE_CSV = Path(config['input'].get('probe_csv')).resolve()

# Genomes
SPECIES_DICT: dict = config.get('species', {})

if not USE_SPECIES_DICT:
    SPECIES: list = [(f.name.split('.fa')[0]).strip() for f in PATH_DICT['SPECIES_DB'].iterdir() if f.suffix == '.fa']
else:
    SPECIES: list = SPECIES_DICT.keys()
