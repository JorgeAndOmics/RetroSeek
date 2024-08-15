import os
import io
from io import StringIO
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import db_utils
import pprint
import cloudpickle
import pickle

import defaults
from utils import pickler, unpickler
from object_class import Object
from colored_logging import colored_logging

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio.Blast import Record

import logging, coloredlogs

import sys

# Your query sequence in FASTA format
query_sequence = ("GAATTGGAGAAGGCACATATTGTGAGGCCAACACACAGTCCCTTTAACTCCCCAGTATGGCCTGTCAAGAAGCCAGATGGGACCTGGT"
                  "GAATGACTGTGGACTATAGGGAGTTAAATAAGGTGACACCACCCTTGCACGCTGCAGTGCCTTCAATACATGATTTAATGGATCATCTGAC"
                  "TGTCCGCCTGGGGACATATCACTATGTAGTGGACTTGGCCAATGCTTTTTTCTCCATTGACATTGCACCCGAGTATCAAGAGCAGTTTGCTTTTA"
                  "CTTGGGATGGGCGGCAATGGACTTTTCAAGTCCTTCCGCAAGGATACTTGCACAGCCCCACTATCTGCCATGGGCTCGTGGCCCAGGATTTGGCACAG"
                  "TGGGATCGCCCATCCTCTGTGGCCTTGTTTCATTATGTTGATGATATTCTATTAACATCTGATTCTCTTTCTGATTTAGAGCAAGCAGCTCCCTCTC"
                  "TTCTCTGCCACCTGAAGTCACGTGGCTGGGCAGTGAATGAGGAAAAGGTCCAAGGCCCTGGCTTATCCGTCAAGCTTTTGGGTGTTGTGTGGTCGGGT"
                  "AAGACAAAGGTTATACCTGAGGCAATTATTGATAAGATACAAGCCTTTCCCCGGCCGACCAAGGTCTCCCAACTGCAGACATATTTGGGTCTGCTAG"
                  "GATATTGGCGGGCGTTTGTGCCCCATTTAGCACAAATGGCAAGGCCCTTGTACAATATGATAAAA")

# Perform RPS-BLAST using qblast with appropriate parameters
result_handle = NCBIWWW.qblast(
    program="blastx",
    database="cdd",
    sequence=query_sequence,
    entrez_query="(\"Reverse Position Specific BLAST\")",
    expect=0.01,
    hitlist_size=50
)

# Save the results to an XML file
with open("rpsblast_results.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
    blast_records = NCBIXML.parse(out_handle)


    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            if not alignment:
                print("No hits found")
                continue
            for hsp in alignment.hsps:
                print("****Alignment****")
                print(f"sequence: {alignment.title}")
                print(f"length: {alignment.length}")
                print(f"e value: {hsp.expect}")
                print(hsp.query[0:75] + "...")
                print(hsp.match[0:75] + "...")
                print(hsp.sbjct[0:75] + "...")