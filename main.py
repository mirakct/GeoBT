from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML, Record


import tkinter as tk

Entrez.email = "karim.salama@campus.tu-berlin.de"


#a function to get the protein sequence from the accession number
def get_protein_sequence(accession_number):
    handle = Entrez.efetch(db="protein", id=accession_number, rettype="fasta", retmode="text")
    record = handle.read()
    handle.close()
    return record

#a function to blast the protein accession number against the nr database


def blast_protein_sequence(protein_sequence):
    result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence,hitlist_size=1)
    blast_record = NCBIXML.read(result_handle)
    return blast_record

# a function using the parseblasttable to parse the blast record and get the alignment

def get_alignment(blast_record):
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print('****Alignment****')
            print('sequence:', alignment.title)
            print('length:', alignment.length)
            print('e value:', hsp.expect)
            print(hsp.query[0:75] + '...')
            print(hsp.match[0:75] + '...')
            print(hsp.sbjct[0:75] + '...')


#test the functions
protein_sequence = get_protein_sequence("NP_001035601.1")
blast_record = blast_protein_sequence(protein_sequence)

#parse the blast record
alignment = get_alignment(blast_record)      
print(blast_record)

#open a tkinter window displaying the alignment


