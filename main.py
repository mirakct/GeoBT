import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
    

def read_fasta(fasta_file):
    """
    Reads a fasta file and returns a list of SeqRecord objects.
    """
    return list(SeqIO.parse(fasta_file, "fasta", generic_dna))

def read_csv(csv_file):
    """
    Reads a csv file and returns a pandas dataframe.
    """
    return pd.read_csv(csv_file)

def write_fasta(fasta_file, seq_records):
    """
    Writes a list of SeqRecord objects to a fasta file.
    """
    SeqIO.write(seq_records, fasta_file, "fasta")

def write_csv(csv_file, dataframe):
    """
    Writes a pandas dataframe to a csv file.
    """
    dataframe.to_csv(csv_file, index=False)

def get_seq(seq_record):
    """
    Returns the sequence of a SeqRecord object.
    """
    return str(seq_record.seq)

def get_id(seq_record):
    """
    Returns the id of a SeqRecord object.
    """
    return seq_record.id
