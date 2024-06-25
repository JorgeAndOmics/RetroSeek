import defaults
from utils import *
import io

test = unpickler(os.path.join('..', 'data', 'pickles'), 'probe_dict.pkl')
# print(test['AJG42161.1'].get_gff())

handle = test['AJG42161.1'].get_fasta()
with tempfile.NamedTemporaryFile(mode='w', delete=True, dir=defaults.TMP_DIR) as fasta_temp_file:
    fasta_temp_file.write(handle)
    fasta_temp_path = fasta_temp_file.name
    print(fasta_temp_path)