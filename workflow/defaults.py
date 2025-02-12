import os

# BLAST
EXPANSION_SIZE: int = 0
E_VALUE: float = 0.09
ACCESSION_ID_REGEX: str = r'[A-Z]{2,}_?[0-9]+\.[0-9]{1,2}'
PROBE_MIN_LENGTH: dict = {
    'GAG': 200,
    'POL': 400,
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
ROOT_CONFIG_FILE = os.path.join('..', 'data', 'config', 'root_folder.txt')
with open(ROOT_CONFIG_FILE, 'r') as f:
    ROOT = f.readline().strip()
ROOT_DB = os.path.join(ROOT, 'local')
SPECIES_DB = os.path.join(ROOT_DB, 'blast_dbs', 'species')
ACCESSORY_DB = os.path.join(ROOT, 'accessory')


# Directories
TABLE_INPUT_DIR = os.path.join('..', 'data', 'tables', 'input')
TABLE_OUTPUT_DIR = os.path.join('..', 'results', 'tables')
TABLE_OVERLAP_MATRIX_DIR = os.path.join(TABLE_OUTPUT_DIR, 'overlap_matrix')
SPECIES_DIR = os.path.join('..', 'data', 'species')
LOG_DIR = os.path.join('..', 'logs')
PICKLE_DIR = os.path.join('..', 'data', 'pickles')
TMP_DIR = os.path.join('..', 'data', 'tmp')
RESULTS_DIR = os.path.join('..', 'results')
PLOT_DIR = os.path.join(RESULTS_DIR, 'plots')
TRACK_DIR = os.path.join(RESULTS_DIR, 'tracks')
TRACK_ORIGINAL_DIR = os.path.join(TRACK_DIR, 'original')
TRACK_CANDIDATES_DIR = os.path.join(TRACK_DIR, 'candidates')
TRACK_VALIDATED_DIR = os.path.join(TRACK_DIR, 'validated')
LTRHARVEST_DIR = os.path.join(RESULTS_DIR, 'ltrharvest')
LTRDIGEST_DIR = os.path.join(RESULTS_DIR, 'ltrdigest')
SEGMENTED_SPECIES_DIR = os.path.join(RESULTS_DIR, 'tables', 'segmented_species')
HMM_PROFILE_DIR = os.path.join(ROOT, 'accessory', 'hmm_profiles')


# Directory generation
os.makedirs(ROOT, exist_ok=True)
os.makedirs(ROOT_DB, exist_ok=True)
os.makedirs(SPECIES_DIR, exist_ok=True)
os.makedirs(SPECIES_DB, exist_ok=True)
os.makedirs(ACCESSORY_DB, exist_ok=True)
os.makedirs(LTRHARVEST_DIR, exist_ok=True)
os.makedirs(LTRDIGEST_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(TRACK_DIR, exist_ok=True)
os.makedirs(SEGMENTED_SPECIES_DIR, exist_ok=True)
os.makedirs(TRACK_ORIGINAL_DIR, exist_ok=True)
os.makedirs(TRACK_CANDIDATES_DIR, exist_ok=True)
os.makedirs(TRACK_VALIDATED_DIR, exist_ok=True)
os.makedirs(TABLE_INPUT_DIR, exist_ok=True)
os.makedirs(TABLE_OUTPUT_DIR, exist_ok=True)
os.makedirs(TABLE_OVERLAP_MATRIX_DIR, exist_ok=True)
os.makedirs(PICKLE_DIR, exist_ok=True)
os.makedirs(HMM_PROFILE_DIR, exist_ok=True)


# Execution and requests
USE_SPECIES_LIST: bool = False
MAX_RETRIEVAL_ATTEMPTS: int = 3
MAX_EXECUTION_ATTEMPTS_PER_SECOND: int = 10
MIN_EXECUTION_INTERVAL: int = 1  # seconds
MAX_THREADPOOL_WORKERS: int = None  # In my laptop, 7 is the maximum number of workers that can be used
GENBANK_RETRIEVAL: bool = True
ENTREZ_EMAIL: str = 'jgonzlez@tcd.ie'
NCBI_API_TOKEN: str = 'faa9e17bb461e82963f079c167ec5c7aac08'


# Displays
DISPLAY_REQUESTS_WARNING: bool = False
DISPLAY_OPERATION_INFO: bool = False


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


# Genomes
SPECIES_FILE = os.path.join(SPECIES_DIR, 'species.txt')

if not USE_SPECIES_LIST:
    SPECIES: list = [(f.split('.fa')[0]).strip() for f in os.listdir(SPECIES_DB) if f.endswith('.fa')]
else:
    SPECIES: list = [line.strip() for line in open(SPECIES_FILE, 'r')]

VIRUS: list = ['NCBI_Virus']

# [
#         'Desmodus_rotundus',
#         'Miniopterus_schreibersii',
#         'Tadarida_brasiliensis',
#         'Antrozous_pallidus',
#         'Molossus_molossus',
#         'Artibeus_lituratus',
#         'Eptesicus_fuscus',
#         'Myotis_myotis',
#         'Eptesicus_nilssonii',
#         'Pipistrellus_kuhlii',
#         'Rhinolophus_ferrumequinum',
#         'Saccopteryx_bilineata',
#         'Vespertilio_murinus',
#         'Plecotus_auritus',
#         'Rhinolophus_hipposideros',
#         'Phyllostomus_discolor',
#         'Myotis_daubentonii',
#         'Myotis_mystacinus',
#         'Corynorhinus_townsendii',
#         'Hipposideros_larvatus',
#         'Rhynchonycteris_naso',
#         'Saccopteryx_leptura',
#         'Molossus_alvarezi',
#         'Glossophaga_mutica',
#         'Molossus_nigricans'
#     ]