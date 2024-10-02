# TODO: Database generation and sorting
from collections import defaultdict
import os

import seq_utils
import utils
import defaults
from object_class import Object

from colored_logging import colored_logging
import logging


def reciprocal_blast(object_dict_a: dict,
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
            :param object_dict_a: The dictionary containing the A (posterior) object pairs.
            :param object_dict_b: The dictionary containing the B (prior) object pairs.
            :param command_a: The BLAST command to run from A to B (e.g. "blastx").
            :param command_b: The BLAST command to run from B to A (e.g. "tblastn").
            :param input_database_a: The path to the dictionary a (e.g. Virus) database.
            :param input_database_b: The path to the dictionary b (e.g. Probes) database.

        Returns:
            :returns: A list of dictionaries containing the reciprocal hits as key value dictionaries.
    """

    # Let's assume databases are already created. [blast_threadpool_executor] already generates
    # fasta files for the input.
    def run_blast_looper(object_dict: dict,
                         command: str,
                         input_database_path) -> dict[list]:
        """
        Run BLAST for a dictionary of objects and return the results in a dictionary.

            Parameters
            ----------
                :param object_dict: The dictionary containing the object pairs.
                :param command: The BLAST command to run (e.g. "blastx").
                :param input_database_path: The path to the input database.

            Returns
            -------
                :returns: A dictionary containing the BLAST results for each object.

        """
        results = defaultdict(list)
        for key, obj in object_dict.items():
            query: dict = {key: obj}
            blast_threadpool_result: dict = seq_utils.blast_threadpool_executor(object_dict=query,
                                                                                command=command,
                                                                                input_database_path=input_database_path)
            if blast_threadpool_result:
                results[obj].extend(list(blast_threadpool_result.values()))
                # The values are lists of BLAST objects (confirmed)!!
            else:
                logging.warning(f'No hits found for {key}.')

        for key, value in results.items():
            logging.info(f'Value type: {type(value)}')
            logging.info(f'Value: {value}')
            for obj in value:
                logging.info(f'Object type: {type(obj)}')
                logging.info(f'Object: {obj}')

        return results

    # TODO: EXPERIMENTAL: Filters a range of returns based on a minimum threshold for bitscore, identity, and E-value.
    def get_range_blast_hit(blast_results: dict, bitscore: int, identity: float, e_value: float) -> dict:
        return {key: [value for value in values if value.HSP.bits >= bitscore and
                      (value.HSP.identities / value.HSP.align_length) >= identity and
                      value.HSP.expect <= e_value] for key, values in blast_results.items()}

    def get_max_blast_hit(blast_results: dict) -> dict:
        """
        Get the top hit for each object in the BLAST results.

            Parameters
            ----------
                :param blast_results: The dictionary containing the BLAST results.

            Returns
            -------
                :returns: A dictionary containing the top hits for each object in the format Object: Object.
        """
        logging.debug('Starting max BLAST hit selection')

        # Get all hits with the maximum bitscore
        selected_hits = {key: [obj for obj in value if float(obj.HSP.bits) == max(float(o.HSP.bits) for o in value)]
                         for key, value in blast_results.items()}

        print(f'Selected hits after bitscore: {selected_hits}')

        # Refine the selection if there are multiple hits with the same bitscore
        for key, value in selected_hits.items():
            if len(value) > 1:
                logging.warning(f'Multiple hits found for {key} with the same bitscore.')

                # Get all hits with the maximum identity score
                max_identity = max(float(obj.HSP.identities) / float(obj.HSP.align_length) for obj in value)
                selected_hits[key] = [obj for obj in value if
                                      (float(obj.HSP.identities) / float(obj.HSP.align_length)) == max_identity]

        # Refine the selection further if there are still multiple hits with the same identity score
        for key, value in selected_hits.items():
            if len(value) > 1:
                logging.warning(f'Multiple hits found for {key} with the same identity score.')

                # Get all hits with the minimum e-value
                min_e_value = min(float(obj.HSP.expect) for obj in value)
                selected_hits[key] = [obj for obj in value if float(obj.HSP.expect) == min_e_value]

        # Finally, select the first hit if there are still multiple hits
        for key, value in selected_hits.items():
            if len(value) > 1:
                logging.warning(f'Multiple hits found for {key} with the same e-value.')

            selected_hits[key] = value[0]  # Select always the first hit at this point to ensure the value is not a list

        return selected_hits

    def hit_description(instance) -> str:
        """
        Get the hit description from the instance from the hit_def attribute.

            Parameters
            ----------
                :param instance: The instance to get the hit description from.

            Returns
            -------
                :returns: The hit description.

        """
        if isinstance(instance, Object):
            return f'{instance.alignment.hit_def.split("(")[0].strip()}'

    # TODO: EXPERIMENTAL: Reciprocal hit finder for a range of returns based on a minimum threshold for bitscore, identity, and E-value.
    def reciprocal_hit_range(route_a_max: dict, route_b_max: dict) -> list[dict]:
        reciprocal_hits = {}
        for x, y in route_a_max.items():  # x: Virus Object w/f; y: Probe Object List
            for i, j in route_b_max.items():  # i: Probe Object w/f; j: Virus Object List
                for y_obj in y:
                    for j_obj in j:
                        if x.species == y.species:
                            if hit_description(x) == hit_description(j_obj) and hit_description(y) == hit_description(
                                    y_obj):
                                logging.info(
                                    f'Walking graph: {x.identifier} | {j_obj.identifier} and '
                                    f'{hit_description(y)} | {hit_description(y_obj)}: Match found')
                                reciprocal_hits |= {x: y}

                            else:
                                logging.debug(
                                    f'Walking graph: {x.identifier} | {j_obj.identifier} and '
                                    f'{hit_description(y)} | {hit_description(y_obj)}')

    # Can't compare instances through identifiers, due to them being randomly generated at BLAST parsing.
    def reciprocal_hit_finder(route_a_max: dict, route_b_max: dict) -> dict:
        """
        Find reciprocal hits between two dictionaries of max BLAST hits.

            Parameters
            ----------
                :param route_a_max: The dictionary containing the max hits for route A.
                :param route_b_max: The dictionary containing the max hits for route B.

            Returns
            -------
                :returns: A list of dictionaries containing the reciprocal hits as key value dictionaries.

        """
        reciprocal_hits: dict = {}
        for x, y in route_a_max.items():  # x: A Object w/f; y: B Object
            for i, j in route_b_max.items():  # i: B Object w/f; j: A Object
                if (x.species == y.species and
                        (hit_description(x) == hit_description(j) and
                         hit_description(y) == hit_description(i))):
                    reciprocal_hits |= {x: i}
                    logging.info(
                        f'Walking graph: {x.identifier} -> {j.identifier} | {hit_description(y)} -> {hit_description(i)} | '
                        f'{hit_description(x)} -> {hit_description(j)}: Match found')
                    # TODO: LTRHARVEST IS ADDING EXTRA DATA TO THE HIT_DEF. THIS IS WHY THE IDENTIFIERS ARE NOT
                    #  MATCHING. HAVE TO UNDERSTAND THE LOGIC AND LOOK FOR BETTER IDENTIFICATION METHODS
                    #  (PERHAPS IN [HIT_DESCRIPTION?]
                    reciprocal_hits |= {x: y}

                else:
                    logging.debug(
                        f'Walking graph: {x.identifier} -> {j.identifier} | {hit_description(y)} -> {hit_description(i)} | '
                        f'{hit_description(x)} -> {hit_description(j)}')


        return reciprocal_hits

    # Run BLAST from A to B
    # Dictionary: Key=Virus object, Value=List of Probe BLAST objects
    route_a: dict = run_blast_looper(object_dict=object_dict_a,
                                     command=command_a,
                                     input_database_path=input_database_b)

    # Run BLAST from B to A
    # Dictionary: Key=Probe object, Value=List of Virus BLAST objects
    route_b: dict = run_blast_looper(object_dict=object_dict_b,
                                     command=command_b,
                                     input_database_path=input_database_a)

    # Get max BLAST hits
    # Dictionary: Key=Virus object, Value=Max Probe BLAST object
    route_a_max: dict = get_max_blast_hit(blast_results=route_a)
    # Dictionary: Key=Probe object, Value=Max Virus BLAST object
    route_b_max: dict = get_max_blast_hit(blast_results=route_b)

    # Find reciprocal hits
    reciprocal_hits: dict = reciprocal_hit_finder(route_a_max=route_a_max,
                                                  route_b_max=route_b_max)

    return reciprocal_hits


if __name__ == '__main__':
    colored_logging(log_file_name='rec_blaster.txt')

    a_dict: dict = utils.unpickler(input_directory_path=defaults.PICKLE_DIR,
                                   input_file_name='full_genome_blast.pkl')

    b_dict: dict = utils.unpickler(input_directory_path=defaults.PICKLE_DIR,
                                   input_file_name='ltr_fasta_blast.pkl')

    reciprocal_hits: list[dict] = reciprocal_blast(object_dict_a=a_dict,
                                                   object_dict_b=b_dict,
                                                   command_a='blastn',
                                                   command_b='blastn',
                                                   input_database_a=os.path.join(defaults.A_END_REC_DB,
                                                                                 'a', 'a'),
                                                   input_database_b=os.path.join(defaults.B_END_REC_DB,
                                                                                 'b', 'b'))

    utils.pickler(data=reciprocal_hits,
                  output_directory_path=defaults.PICKLE_DIR,
                  output_file_name='reciprocal_hits.pkl')

    logging.info(f'Successfully performed reciprocal BLAST on {len(reciprocal_hits)} sequences.')
