import yaml
import os

CONFIG_FILE = os.path.abspath(os.path.join('..', 'data', 'config', 'config.yaml'))
with open(CONFIG_FILE, 'r') as f:
    config = yaml.safe_load(f)

# BLAST
EXPANSION_SIZE: int = 0
E_VALUE = config['blast']['e_value']
ACCESSION_ID_REGEX = config['blast']['accession_id_regex']
PROBE_MIN_LENGTH = config['blast']['probe_min_length']

# Logging
LEVEL_STYLES = config['logging']['level_styles']
FIELD_STYLES = config['logging']['field_styles']

# Databases
ROOT = os.path.abspath(config['root'].get('db_root_folder', os.path.join('..')))
ROOT_DB = os.path.abspath(os.path.join(ROOT, 'local'))
SPECIES_DB = os.path.abspath(os.path.join(ROOT_DB, 'blast_dbs', 'species'))
ACCESSORY_DB = os.path.abspath(os.path.join(ROOT, 'accessory'))

# Directories
DATA_DIR = os.path.abspath(os.path.join(config['root'].get('data_root_folder', os.path.join(ROOT, 'data'))))
RESULTS_DIR = os.path.abspath(os.path.join(config['root'].get('results_root_folder', os.path.join(ROOT, 'results'))))
LOG_DIR = os.path.abspath(os.path.join(config['root'].get('log_root_folder', os.path.join(ROOT, 'logs'))))
WORKFLOW_DIR = os.path.abspath(os.path.join(config['root'].get('workflow_root_folder', os.path.join(ROOT, 'workflow'))))

DOWNLOAD_LOG = os.path.abspath(os.path.join(LOG_DIR, 'download_log.log'))

HMM_PROFILE_DIR = os.path.abspath(os.path.join(ACCESSORY_DB, 'hmm_profiles'))

SPECIES_DIR = os.path.abspath(os.path.join(DATA_DIR, 'species'))
TABLE_INPUT_DIR = os.path.abspath(os.path.join(DATA_DIR, 'tables'))
PICKLE_DIR = os.path.abspath(os.path.join(DATA_DIR, 'pickles'))
TMP_DIR = os.path.abspath(os.path.join(DATA_DIR, 'tmp'))

TBLASTN_PICKLE_DIR = os.path.abspath(os.path.join(PICKLE_DIR, 'tblastn'))

TABLE_OUTPUT_DIR = os.path.abspath(os.path.join(RESULTS_DIR, 'tables'))
PLOT_DIR = os.path.abspath(os.path.join(RESULTS_DIR, 'plots'))
TRACK_DIR = os.path.abspath(os.path.join(RESULTS_DIR, 'tracks'))
LTRHARVEST_DIR = os.path.abspath(os.path.join(RESULTS_DIR, 'ltrharvest'))
LTRDIGEST_DIR = os.path.abspath(os.path.join(RESULTS_DIR, 'ltrdigest'))

CIRCLE_PLOT_DIR = os.path.abspath(os.path.join(PLOT_DIR, 'circle_plots'))
HOTSPOT_PDF_DIR = os.path.abspath(os.path.join(PLOT_DIR, 'hotspot_pdfs'))

TRACK_ORIGINAL_DIR = os.path.abspath(os.path.join(TRACK_DIR, 'original'))
TRACK_CANDIDATES_DIR = os.path.abspath(os.path.join(TRACK_DIR, 'candidates'))
TRACK_VALIDATED_DIR = os.path.abspath(os.path.join(TRACK_DIR, 'validated'))
TRACK_HOTSPOTS_DIR = os.path.abspath(os.path.join(TRACK_DIR, 'hotspots'))

TABLE_HOTSPOT_DIR = os.path.abspath(os.path.join(TABLE_OUTPUT_DIR, 'hotspots'))
TABLE_OVERLAP_MATRIX_DIR = os.path.abspath(os.path.join(TABLE_OUTPUT_DIR, 'overlap_matrix'))
SEGMENTED_SPECIES_DIR = os.path.abspath(os.path.join(TABLE_OUTPUT_DIR, 'segmented_species'))
PLOT_DATAFRAMES_DIR = os.path.abspath(os.path.join(TABLE_OUTPUT_DIR, 'plot_dataframes'))

# Directory generation
os.makedirs(ROOT, exist_ok=True)
os.makedirs(ROOT_DB, exist_ok=True)
os.makedirs(SPECIES_DIR, exist_ok=True)
os.makedirs(SPECIES_DB, exist_ok=True)
os.makedirs(ACCESSORY_DB, exist_ok=True)
os.makedirs(LTRHARVEST_DIR, exist_ok=True)
os.makedirs(LTRDIGEST_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)
os.makedirs(CIRCLE_PLOT_DIR, exist_ok=True)
os.makedirs(HOTSPOT_PDF_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(TBLASTN_PICKLE_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(WORKFLOW_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(TRACK_DIR, exist_ok=True)
os.makedirs(SEGMENTED_SPECIES_DIR, exist_ok=True)
os.makedirs(TRACK_ORIGINAL_DIR, exist_ok=True)
os.makedirs(TRACK_CANDIDATES_DIR, exist_ok=True)
os.makedirs(TRACK_VALIDATED_DIR, exist_ok=True)
os.makedirs(TRACK_HOTSPOTS_DIR, exist_ok=True)
os.makedirs(TABLE_INPUT_DIR, exist_ok=True)
os.makedirs(TABLE_OUTPUT_DIR, exist_ok=True)
os.makedirs(TABLE_OVERLAP_MATRIX_DIR, exist_ok=True)
os.makedirs(TABLE_HOTSPOT_DIR, exist_ok=True)
os.makedirs(PICKLE_DIR, exist_ok=True)
os.makedirs(PLOT_DATAFRAMES_DIR, exist_ok=True)
os.makedirs(HMM_PROFILE_DIR, exist_ok=True)

# Execution and requests
USE_SPECIES_LIST = config['execution']['use_species_list']
MAX_RETRIEVAL_ATTEMPTS = config['execution']['max_retrieval_attempts']
MAX_THREADPOOL_WORKERS = config['execution']['max_threadpool_workers']
ENTREZ_EMAIL = config['execution']['entrez_email']

# Display
DISPLAY_REQUESTS_WARNING: bool = config['display']['display_requests_warning']
DISPLAY_OPERATION_INFO: bool = config['display']['display_operation_info']

# CSV
CSV_DELIMITER: str = config['tables']['delimiter']
CSV_ATTRIBUTES: list[str] = config['tables']['attributes']

# Genomes
SPECIES_FILE: str = config['input']['species_file']

if not USE_SPECIES_LIST:
    SPECIES: list = [(f.split('.fa')[0]).strip() for f in os.listdir(SPECIES_DB) if f.endswith('.fa')]
else:
    SPECIES: list = [line.strip() for line in open(SPECIES_FILE, 'r')]
