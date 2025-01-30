import os

# GENOME DOWNLOAD
ASSEMBLY_PRIORITY = {  # Define an ordering for assembly levels from highest to lowest
    "Complete Genome": 1,
    "Chromosome": 2,
    "Scaffold": 3,
    "Contig": 4
}

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
ROOT = os.path.join('/mnt/', 'v', 'databases')
ROOT_DB = os.path.join(ROOT, 'local')
SPECIES_DB = os.path.join(ROOT_DB, 'blast_dbs', 'species')
VIRUS_DB = os.path.join(ROOT_DB, 'blast_dbs', 'virus')
A_END_REC_DB = os.path.join(ROOT_DB, 'rec_dbs', 'a_point_rec')
B_END_REC_DB = os.path.join(ROOT_DB, 'rec_dbs', 'b_point_rec')
LTR_DB = os.path.join(ROOT_DB, 'ltr_dbs')

# Directories
TABLE_INPUT_DIR = os.path.join('..', 'data', 'tables', 'input')
TABLE_OUTPUT_DIR = os.path.join('..', 'results', 'tables')
TABLE_OVERLAP_MATRIX_DIR = os.path.join(TABLE_OUTPUT_DIR, 'overlap_matrix')
LOG_DIR = os.path.join('..', 'logs')
PICKLE_DIR = os.path.join('..', 'data', 'pickles')
TMP_DIR = os.path.join('..', 'data', 'tmp')
RESULTS_DIR = os.path.join('..', 'results')
TRACK_DIR = os.path.join(RESULTS_DIR, 'tracks')
TRACK_CANDIDATES_DIR = os.path.join(TRACK_DIR, 'candidates')
TRACK_VALIDATED_DIR = os.path.join(TRACK_DIR, 'validated')
LTRHARVEST_DIR = os.path.join(RESULTS_DIR, 'ltrharvest')
LTRDIGEST_DIR = os.path.join(RESULTS_DIR, 'ltrdigest')
SEGMENTED_SPECIES_DIR = os.path.join(RESULTS_DIR, 'tables', 'segmented_species')
HMM_PROFILE_DIR = os.path.join(ROOT, 'accessory', 'hmm_profiles')

# Execution and requests
MAX_RETRIEVAL_ATTEMPTS: int = 3
MAX_EXECUTION_ATTEMPTS_PER_SECOND: int = 10
MIN_EXECUTION_INTERVAL: int = 1  # seconds
MAX_THREADPOOL_WORKERS: int = None  # In laptop, 7 is the maximum number of workers that can be used
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
SPECIES: list = [
                   'Desmodus_rotundus',
                   'Miniopterus_schreibersii',
                   'Tadarida_brasiliensis',
                   'Antrozous_pallidus',
                   'Molossus_molossus',
                   'Artibeus_lituratus',
                   'Eptesicus_fuscus',
                     'Myotis_myotis',   #faulty
                     'Eptesicus_nilssonii',
                     'Pipistrellus_kuhlii',   #faulty
                     'Rhinolophus_ferrumequinum',
                     'Saccopteryx_bilineata',
                     'Vespertilio_murinus',
                     'Plecotus_auritus',
                     'Rhinolophus_hipposideros',
                     'Phyllostomus_discolor',
                     'Myotis_daubentonii',
                     'Myotis_mystacinus',
                     'Corynorhinus_townsendii',
                     'Hipposideros_larvatus',
                     'Rhynchonycteris_naso',
                     'Saccopteryx_leptura',
                     'Molossus_alvarezi',
                     'Glossophaga_mutica',
                     'Molossus_nigricans'   #faulty
                 ]

VIRUS: list = ['NCBI_Virus']

