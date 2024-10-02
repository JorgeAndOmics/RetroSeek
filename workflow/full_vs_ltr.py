import seq_utils
import utils
import defaults

from colored_logging import colored_logging
import logging




if __name__ == '__main__':
    colored_logging(log_file_name='full_vs_ltr.txt')

    full_genome_dict: dict = utils.unpickler(input_directory_path=defaults.PICKLE_DIR,
                                             input_file_name='full_genome_blast.pkl')

    # Divide full_genome_dict into different dictionaries based on the species contained in the objects
    # This is necessary because the LTR BLAST database is also divided by species
    species_dict: dict = seq_utils.species_divider(object_dict=full_genome_dict)

    result_dict: dict = {}
    # Process each species as a separate dictionary
    for species, species_objects in species_dict.items():
        if species in defaults.SPECIES:
            logging.debug(f'Performing BLASTn and retrieval on {species}')
            # Perform tBLASTn and retrieval on the species objects
            blastn2gb_results: dict = seq_utils.blast_retriever(object_dict=species_objects,
                                                                command='blastn',
                                                                genome=[species],  # Pass the single species in loop as a
                                                                # list of one element
                                                                online_database='nucleotide',
                                                                input_database_path=defaults.LTR_DB,
                                                                multi_threading=True)

            result_dict |= blastn2gb_results

    utils.pickler(
        data=result_dict,
        output_directory_path=defaults.PICKLE_DIR,
        output_file_name='full_vs_ltr.pkl',
    )

    # utils.directory_content_eraser(directory_path=defaults.TMP_DIR)

    logging.info(f'Successfully performed BLASTn and retrieval on {len(result_dict)} sequences.')
