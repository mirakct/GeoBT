# Import important packages
from Bio import Entrez
from Bio import SeqIO   

# Download protein from NCBI 
def download_protein_sequence(accession):                                   # Defining the function
    Entrez.email = 'linda.jungwirth1@gmail.com'                             # Set your email address for Entrez
    handle = Entrez.efetch(db='protein',                                    # Defining the handle
                           id=accession,                                    # Defining the accession number
                           rettype='fasta',                                 # Defining the return type
                           retmode='text')                                  # Defining the return mode
    return handle.read()                                                    # Returning the record

# Download core-MetE protein sequence
accession = 'WP_011309033'                                                  # Defining the accession number
sequence = download_protein_sequence(accession)                             # Defining the sequence
print(sequence)                                                             # Printing the sequence

#Save sequence as a fasta file  #needs the SeqIO module
file_path = "/Users/lindajungwirth/python/sequence.fasta"                   # Defining the file path
with open(file_path, "w") as file:                                          # Opening the file
    file.write(sequence)                                                    # Writing the sequence to the file
print("Saved sequence to", file_path)                                       # Printing the file path