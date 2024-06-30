# TODO: Database generation and sorting
from collections import defaultdict

from seq_utils import *
from utils import *

from colored_logging import colored_logging
import logging

def reciprocal_blast_experimental(object_dict_a: dict,
                                  object_dict_b: dict,
                                  command_a: str,
                                  command_b: str,
                                  input_database_a: str,
                                  input_database_b: str):
    """
    Perform reciprocal BLAST between two object pair dictionaries, provided existing databases from their FASTA files.
    route_a: virus -> probes
    route_b: probes -> virus

        Args:
            object_dict_a (dict): The dictionary containing the Virus object pairs.
            object_dict_b (dict): The dictionary containing the Probes object pairs.
            command_a (str): The BLAST command to run from A to B (e.g. "blastx").
            command_b (str): The BLAST command to run from B to A (e.g. "tblastn").
            input_database_a (str): The path to the dictionary a (e.g. Virus) database.
            input_database_b (str): The path to the dictionary b (e.g. Probes) database.

        Returns:
            list: A list of dictionaries containing the reciprocal hits as key value dictionaries.
    """

    # Let's assume databases are already created. [blast_threadpool_executor] already generates
    # fasta files for the input.
    def run_blast_looper(object_dict, command, input_database_path):
        results = defaultdict(list)
        for key, obj in object_dict.items():
            query = {key: obj}
            results[key].append(blast_threadpool_executor(query, command, input_database_path).values())
        return results

    def get_max_blast_hit(blast_results):
        # Ordered based on HSP bitscore
        return {key: max(value, key=lambda x: x.HSP.bits) for key, value in blast_results.items()}

    def fasta_desc(instance):
        if isinstance(instance, Object):
            return instance.get_fasta('seqrecord').description

    # Can't compare instances through identifiers, due to them being randomly generated at BLAST parsing.
    def reciprocal_hit_finder(route_a_max, route_b_max):
        reciprocal_hits = []
        for x, y in route_a_max.items():  # x: Virus Object; y: Probe Object
            for i, j in route_b_max.items():  # i: Probe Object; j: Virus Object
                if fasta_desc(x) == fasta_desc(j) and fasta_desc(y) == fasta_desc(i):  # Now it compares FASTA headers
                # directly. HSP selection is made at max() bitscore selection.
                    reciprocal_hits.append({x: y})
        return reciprocal_hits

    # Run BLAST from A to B
    # Dictionary: Key=Virus object, Value=List of Probe BLAST objects
    route_a: dict = run_blast_looper(object_dict=object_dict_a,
                                     command=command_a,
                                     input_database_path=input_database_a)

    # Run BLAST from B to A
    # Dictionary: Key=Probe object, Value=List of Virus BLAST objects
    route_b: dict = run_blast_looper(object_dict=object_dict_b,
                                     command=command_b,
                                     input_database_path=input_database_b)

    # Get max BLAST hits
    # Dictionary: Key=Virus object, Value=Max Probe BLAST object
    route_a_max = get_max_blast_hit(blast_results=route_a)
    # Dictionary: Key=Probe object, Value=Max Virus BLAST object
    route_b_max = get_max_blast_hit(blast_results=route_b)

    # Find reciprocal hits
    reciprocal_hits = reciprocal_hit_finder(route_a_max=route_a_max,
                                           route_b_max=route_b_max)

    return reciprocal_hits


if __name__ == '__main__':
    colored_logging(log_file_name='rec_blaster.txt')

    virus_dict = unpickler(input_directory_path=os.path.join('..', 'data', 'pickles'),
                           input_file_name='virus_dict.pkl')

    probe_dict = unpickler(input_directory_path=os.path.join('..', 'data', 'pickles'),
                           input_file_name='probe_dict.pkl')

    reciprocal_hits = reciprocal_blast_experimental(object_dict_a=virus_dict,
                                                    object_dict_b=probe_dict,
                                                    command_a='blastx',
                                                    command_b='tblastn',
                                                    input_database_a='placeholder',
                                                    input_database_b='placeholder')