# TODO: Database generation and sorting
from collections import defaultdict

from workflow import seq_utils
from workflow import utils
from workflow import defaults
from workflow.object_class import Object

from workflow.colored_logging import colored_logging
import logging


def reciprocal_blast_experimental(object_dict_a: dict,
                                  object_dict_b: dict,
                                  command_a: str,
                                  command_b: str,
                                  input_database_a: str,
                                  input_database_b: str) -> list[dict]:
    """
    Perform reciprocal BLAST between two object pair dictionaries, provided existing databases from their FASTA files.
    route_a: virus -> probes
    route_b: probes -> virus

        Parameters
        ------------
            :param object_dict_a: The dictionary containing the Virus object pairs.
            :param object_dict_b: The dictionary containing the Probes object pairs.
            :param command_a: The BLAST command to run from A to B (e.g. "blastx").
            :param command_b: The BLAST command to run from B to A (e.g. "tblastn").
            :param input_database_a: The path to the dictionary a (e.g. Virus) database.
            :param input_database_b: The path to the dictionary b (e.g. Probes) database.

        Returns:
            :returns: A list of dictionaries containing the reciprocal hits as key value dictionaries.
    """

    # Let's assume databases are already created. [blast_threadpool_executor] already generates
    # fasta files for the input.
    def run_blast_looper(object_dict, command, input_database_path) -> dict[list]:
        results = defaultdict(list)
        for key, obj in object_dict.items():
            query: dict = {key: obj}
            results[key].append(seq_utils.blast_threadpool_executor(object_dict=query,
                                                                    command=command,
                                                                    input_database_path=input_database_path).values())
        return results

    # TODO: EXPERIMENTAL: Filters a range of returns based on a minimum threshold for bitscore, identity, and E-value.
    def get_range_blast_hit(blast_results: dict, bitscore: int, identity: float, e_value: float) -> dict:
        return {key: [value for value in values if value.HSP.bits >= bitscore and
                      (value.HSP.identities / value.HSP.align_length) >= identity and
                      value.HSP.expect <= e_value] for key, values in blast_results.items()}

    def get_max_blast_hit(blast_results) -> dict:
        # Ordered first based on HSP bitscore
        selected_hits: dict = {key: max(value, key=lambda x: x.HSP.bits) for key, value in blast_results.items()}
        if len(selected_hits.values()) > 1:
            logging.warning('Multiple hits found with the same bitscore.')
            selected_hits: dict = {key: max(selected_hits, key=lambda x: (x.HSP.identities / x.HSP.align_length)) for
                                   key, value in selected_hits.items()}
            if len(selected_hits.values()) > 1:
                logging.warning('Multiple hits found with the same identity score.')
                selected_hits: dict = {key: min(selected_hits, key=lambda x: x.HSP.expect) for
                                       key, value in blast_results.items()}
                if len(selected_hits.values()) > 1:
                    logging.warning('Multiple hits found with same E-value.')
                    selected_hits: dict = {key: value[0] for key, value in selected_hits.items}

        return selected_hits

    def hit_desc(instance) -> str:
        if isinstance(instance, Object):
            return instance.alignment.hit_def

    # TODO: EXPERIMENTAL: Reciprocal hit finder for a range of returns based on a minimum threshold for bitscore, identity, and E-value.
    def reciprocal_hit_range(route_a_max: dict, route_b_max: dict) -> list[dict]:
        reciprocal_hits = []
        for x, y in route_a_max.items():  # x: Virus Object w/f; y: Probe Object List
            for i, j in route_b_max.items():  # i: Probe Object w/f; j: Virus Object List
                for y_obj in y:
                    for j_obj in j:
                        if hit_desc(x) == hit_desc(j_obj) and hit_desc(y) == hit_desc(y_obj):
                            logging.info(
                                f'Walking graph: {x.identifier} | {j_obj.identifier} and '
                                f'{hit_desc(y)} | {hit_desc(y_obj)}: Match found')
                            reciprocal_hits.append({x: y})

                        else:
                            logging.debug(
                                f'Walking graph: {x.identifier} | {j_obj.identifier} and '
                                f'{hit_desc(y)} | {hit_desc(y_obj)}')

    # Can't compare instances through identifiers, due to them being randomly generated at BLAST parsing.
    def reciprocal_hit_finder(route_a_max: dict, route_b_max: dict) -> list[dict]:
        reciprocal_hits = []
        for x, y in route_a_max.items():  # x: Virus Object w/f; y: Probe Object
            for i, j in route_b_max.items():  # i: Probe Object w/f; j: Virus Object
                if hit_desc(x) == hit_desc(j) and hit_desc(y) == hit_desc(i):  # Now it compares alignment.hit_def
                    # attributes directly (identical to FASTA headers). HSP selection is made at max() metric selection.
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
    route_a_max: dict = get_max_blast_hit(blast_results=route_a)
    # Dictionary: Key=Probe object, Value=Max Virus BLAST object
    route_b_max: dict = get_max_blast_hit(blast_results=route_b)

    # Find reciprocal hits
    reciprocal_hits: list[dict] = reciprocal_hit_finder(route_a_max=route_a_max,
                                                        route_b_max=route_b_max)

    return reciprocal_hits


if __name__ == '__main__':
    colored_logging(log_file_name='rec_blaster.txt')

    virus_dict: dict = utils.unpickler(input_directory_path=defaults.PICKLE_DIR,
                                       input_file_name='virus_dict.pkl')

    probe_dict: dict = utils.unpickler(input_directory_path=defaults.PICKLE_DIR,
                                       input_file_name='probe_dict.pkl')

    reciprocal_hits: list[dict] = reciprocal_blast_experimental(object_dict_a=virus_dict,
                                                                object_dict_b=probe_dict,
                                                                command_a='blastx',
                                                                command_b='tblastn',
                                                                input_database_a='placeholder',
                                                                input_database_b='placeholder')