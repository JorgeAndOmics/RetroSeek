# TODO: COMPLETELY TEMPORARY, FINAL IMPLEMENTATION WILL USE AUTOMATED EXISTING BIOINFORMATICS TOOLS

import csv
import logging
import os
import urllib.request
import gzip
import shutil
import defaults

from Bio import Entrez
from Bio import SeqIO  # Import SeqIO for parsing
from colored_logging import colored_logging
from collections import Counter

def read_tax_ids_from_csv(file_path: os.PathLike, column_name: str):
    """
    Reads a CSV file and extracts a column containing tax-ids.
    """
    tax_ids = []
    with open(file_path, mode='r') as file:
        logging.debug('Reading tax-ids from CSV file')
        reader = csv.DictReader(file)
        tax_ids.extend(row[column_name] for row in reader)
    return [int(tax_id) for tax_id in tax_ids if tax_id.isdigit()]


def get_assembly_summary(assembly_id: str):
    """
    Retrieves the assembly summary from Entrez for a given assembly ID.
    """
    handle = Entrez.esummary(db="assembly", id=assembly_id)
    summary = Entrez.read(handle)
    return summary['DocumentSummarySet']['DocumentSummary'][0]


def get_best_assembly_accession(tax_id: int, assembly_priority: dict = defaults.ASSEMBLY_PRIORITY) -> tuple:
    """
    Query NCBI for the best assembly accession using a taxonomic ID.
    """
    Entrez.email = defaults.ENTREZ_EMAIL
    Entrez.api_key = defaults.NCBI_API_TOKEN
    handle = Entrez.esearch(db="assembly", term=f"txid{tax_id}[Organism:exp]", retmax=5)
    record = Entrez.read(handle)

    best_assembly = None
    best_assembly_level = float('inf')  # Start with the worst possible level

    for assembly_id in record['IdList']:
        assembly_info = get_assembly_summary(assembly_id)

        assembly_level = assembly_info['AssemblyStatus']
        if assembly_level in assembly_priority:
            level_priority = assembly_priority[assembly_level]
            if level_priority < best_assembly_level:
                best_assembly_level = level_priority
                ftp_link = assembly_info['FtpPath_GenBank']  # For GenBank assemblies
                best_assembly = (assembly_info['AssemblyAccession'], ftp_link)

    return best_assembly


def genome_exists(tax_id: int, output_dir: os.PathLike) -> bool:
    """
    Checks if the genome assembly for the given tax_id already exists.
    """
    fasta_output_path = os.path.join(output_dir, f"{tax_id}.fna")
    return os.path.exists(fasta_output_path)


def get_assembly_download_link(tax_id: int) -> tuple:
    """
    Gets the download link and assembly accession for the best assembly for the given tax_id.
    """
    best_assembly = get_best_assembly_accession(tax_id)
    if not best_assembly:
        logging.warning(f"No assembly found for tax-id {tax_id}")
        return None, None
    assembly_accession, ftp_link = best_assembly
    gz_file_link = f"{ftp_link}/{ftp_link.split('/')[-1]}_genomic.fna.gz"  # Compressed .gz file
    return assembly_accession, gz_file_link


def download_file_with_retries(url: str, output_path: os.PathLike, retries: int):
    """
    Downloads a file from the given URL to the output path, with the specified number of retries.
    """
    for attempt in range(retries):
        try:
            logging.debug(f"Attempting to download (Attempt {attempt + 1}) from {url}")
            urllib.request.urlretrieve(url, output_path)
            logging.info(f"Saved file to {output_path}")
            return True
        except Exception as e:
            logging.warning(f"Failed to download from {url} on attempt {attempt + 1}: {e}")
            if attempt == retries - 1:
                logging.error(f"Failed to download from {url} after {retries} attempts. Skipping.")
                return False


def uncompress_gz_file(input_path: os.PathLike, output_path: os.PathLike, retries: int):
    """
    Uncompresses a .gz file to the specified output path, with the specified number of retries.
    """
    for attempt in range(retries):
        try:
            with gzip.open(input_path, 'rb') as f_in:
                with open(output_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            logging.info(f"Uncompressed file to {output_path}")
            return True
        except Exception as e:
            logging.warning(f"Failed to uncompress {input_path} on attempt {attempt + 1}: {e}")
            if attempt == retries - 1:
                logging.error(f"Failed to uncompress {input_path} after {retries} attempts. Skipping.")
                return False


def download_and_uncompress_genome_assembly(tax_id: int, output_dir: os.PathLike, retries: int = 3):
    """
    Downloads the genome assembly as a compressed .gz file for a given tax-id,
    uncompresses it, verifies by parsing with SeqIO, and renames the resulting file as {tax_id}.fna.
    A retry mechanism and a check for already downloaded files are included.
    """
    fasta_output_path = os.path.join(output_dir, f"{tax_id}.fna")
    gz_output_path = os.path.join(output_dir, f"{tax_id}.fna.gz")  # Save .gz file

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if genome_exists(tax_id, output_dir):
        logging.info(f"Genome assembly for tax-id {tax_id} already exists. Skipping download.")
        return

    assembly_accession, gz_file_link = get_assembly_download_link(tax_id)
    if not assembly_accession or not gz_file_link:
        return  # Message already logged in get_assembly_download_link

    for attempt in range(retries):
        logging.debug(f"Starting attempt {attempt + 1} for tax-id {tax_id}")
        # Download the .gz file
        if not download_file_with_retries(gz_file_link, gz_output_path, retries=1):
            continue  # Skip to next attempt if download failed

        # Uncompress the .gz file
        if not uncompress_gz_file(gz_output_path, fasta_output_path, retries=1):
            continue  # Skip to next attempt if uncompress failed

        # Verify the uncompressed file by attempting to parse it
        try:
            logging.debug(f"Attempting to parse {fasta_output_path}")
            with open(fasta_output_path, 'r') as handle:
                records = list(SeqIO.parse(handle, 'fasta'))
            if records:
                logging.info(f"Successfully parsed {fasta_output_path}")
                # Optionally remove the .gz file after uncompressing and verifying
                # os.remove(gz_output_path)
                return  # Successfully downloaded, uncompressed, and verified
            else:
                logging.warning(f"No records found in {fasta_output_path}")
        except Exception as e:
            logging.warning(f"Failed to parse {fasta_output_path} on attempt {attempt + 1}: {e}")

        logging.warning(f"Redownloading and uncompressing genome for tax-id {tax_id} (Attempt {attempt + 1})")
        # Remove the possibly corrupted files before retrying
        if os.path.exists(gz_output_path):
            os.remove(gz_output_path)
        if os.path.exists(fasta_output_path):
            os.remove(fasta_output_path)

    logging.error(f"Failed to download and verify genome for tax-id {tax_id} after {retries} attempts. Skipping.")


def download_genomes_for_tax_ids(tax_id_list: list,
                                 output_dir: os.PathLike,
                                 max_attempts: int = defaults.MAX_RETRIEVAL_ATTEMPTS):
    """
    Downloads and uncompresses genome assemblies for all tax-ids in the provided list.
    Includes a retry mechanism and skips files that have already been downloaded.
    """
    total_genomes = len(tax_id_list)
    counter = Counter({'remaining': total_genomes})

    for idx, tax_id in enumerate(tax_id_list, start=1):
        logging.info(f"Downloading genome for tax-id {tax_id} ({idx}/{total_genomes}). {counter['remaining']} genomes remaining.")
        download_and_uncompress_genome_assembly(tax_id,
                                                output_dir,
                                                max_attempts)
        counter['remaining'] -= 1


if __name__ == '__main__':
    colored_logging(log_file_name='genome_downloader.txt')

    # Read tax-ids from CSV, download genomes for all of them
    csv_file_path = '/mnt/f/Descargas/chiroptera_genomes.csv'
    output_path = '/mnt/f/Volumen/genomes_goat_revisited'
    column_name = "taxon_id"  # Column name containing tax-ids in the CSV file

    tax_ids = read_tax_ids_from_csv(csv_file_path,
                                    column_name)

    # Download genomes for the tax-ids from the CSV with retry mechanism
    download_genomes_for_tax_ids(tax_ids,
                                 output_path,
                                 max_attempts=500)
