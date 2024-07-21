from Bio import Entrez, SeqIO
import subprocess
import os

import defaults
import utils

import logging

def get_species_name_from_file(taxid_file: str,
                               _entrez_email: str = defaults.ENTREZ_EMAIL,
                               _entrez_api_token: str = defaults.NCBI_API_TOKEN) -> str:
    """
    Given a taxid number, return the scientific name of the species in the format Genus_species. Used
    in [taxlist2name]. If the file name is not a valid taxid, return the same file name.
    If the queried taxid is not found, return None.

        Parameters
        ----------
            :param taxid_file: Taxonomic ID-named file.
            :param _entrez_email: Email address for Entrez.
            :param _entrez_api_token: API token for Ent.

        Returns
        -------
            :returns: Scientific name in the format Genus_species

    """
    Entrez.email = _entrez_email
    Entrez.api_key = _entrez_api_token
    taxid_file_query = taxid_file.split('.')[0]
    try:
        taxid_file_query = int(taxid_file_query)
        handle = Entrez.efetch(db="taxonomy", id=str(taxid_file_query), retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        if not records:
            logging.warning(f'Taxid {taxid_file} not found.')
            return None

        scientific_name = records[0]['ScientificName']
        return scientific_name.replace(" ", "_")

    except ValueError:
        logging.warning(f'Invalid taxid {taxid_file_query}. Generating folder with the same name.')
        return taxid_file_query



def taxlist2name(taxid_list: list) -> list:
    """
    Given a list of taxids, return the scientific names of the species in the format Genus_species.

        Parameters
        ----------
            :param taxid_list: List of taxonomic IDs.

        Returns
        -------
            :returns: List of scientific names in the format Genus_species.

    """
    return [get_species_name_from_file(taxid) for taxid in taxid_list]


def blast_db_generator(input_file_path,
                       output_directory_path,
                       db_name: str,
                       db_type: str) -> None:
    """
    Generates the BLAST databases from the FASTA files in the specified directory.

        Parameters
        ----------
            :param input_file_path: The input FASTA file.
            :param output_directory_path: The output directory path for the BLAST database.
            :param db_name: The name of the database.
            :param db_type: The type of database to generate.


    """
    makeblastdb_command = ['makeblastdb',
                           '-in', input_file_path,
                           '-dbtype', db_type,
                           '-out', os.path.join(output_directory_path, db_name)]

    subprocess.run(makeblastdb_command)
    logging.info(f'Generated BLAST database for {db_name}.')

def directory_db_generator(file_list: list,
                           input_db,
                           db_type: str) -> list:
    """
    Generates the BLAST databases from the FASTA files in the specified directory.

        Parameters
        ----------
            :param file_list: The list of files in the input directory (tax-id named).
            :param output_list: The (defaults) variable containing the genomes for analysis
            :param input_db: The directory path for the BLAST database (species, virus...).
            :param db_type: The type of database to generate (nucl, prot).

        Returns
        -------
            :return:

    """
    genomes = []
    for file in file_list:  # The original FASTA file paths, with the taxid as the name
        if not os.path.isdir(file):
            file_sn = get_species_name_from_file(file.split('.')[0])
            directory = utils.directory_generator(input_db, file_sn)
            blast_db_generator(input_file_path=os.path.join(input_db, file),
                               output_directory_path=directory,
                               db_name=file_sn,
                               db_type=db_type)

            genomes.extend(file_sn)

    return genomes

def existing_db_list_cleaner(file_list: list, directory_to_check: str) -> list:
    """
    Cleans the file list if the directory with the database already exists, identified by the presence of
    the BLAST database files.

        Parameters
        ----------
            :param file_list: The list of files (tax-id named) in the input directory.
            :param directory_to_check: The directory to check for the database (species, virus...).

        Returns
        -------
            :returns: The cleaned list of files.

    """
    blast_extensions = {'.nhr', '.nin', '.nsq'}
    for dir in os.listdir(directory_to_check):
        full_dir_path = os.path.join(directory_to_check, dir)
        if os.path.isdir(full_dir_path):
            db_exists = any(
                any(file.endswith(ext) for ext in blast_extensions)
                for file in os.listdir(full_dir_path)
            )
            if db_exists:
                logging.info(f'Database {dir} already exists. Cleaning file list.')
                if dir in file_list:
                    file_list.remove(dir)
    return file_list


def objdict2fasta(object_dict: dict,
                  output_directory_path: os.PathLike,
                  output_file_name: str) -> None:
    """
    Converts the object dictionary to a single concatenated FASTA file in the specified directory.

        Parameters
        ----------
            :param object_dict: The dictionary containing the object pairs.
            :param output_directory_path: The output directory path for the FASTA file.
            :param output_filename: The name of the output concatenated FASTA file.
    """
    output_file_path = os.path.join(output_directory_path, output_file_name)

    with open(output_file_path, 'w') as concatenated_fasta:
        for key, obj in object_dict.items():
            SeqIO.write(obj.get_fasta('seqrecord'), concatenated_fasta, 'fasta')
            logging.info(f'Extracted FASTA from {key} and appended to {output_file_name}.')
