# Import important packages
from Bio import Entrez                                                              # Import Entrez
from Bio import SeqIO                                                               # Import SeqIO

# Provide your email address to NCBI
Entrez.email = "jungwirth@campus.tu-berlin.de"

#Get accession number
def get_accession_number():                                                             # Defining the function
    search_term = "Dehalococcoides mccartyi[Organism] AND methionine synthase[Protein]" # Defining the search term
    handle = Entrez.esearch(db="protein", term=search_term, retmax=1)                   # Defining the handle
    record = Entrez.read(handle)                                                        # Reading the handle
    handle.close()                                                                      # Closing the handle

    if int(record["Count"]) == 0:                                                       # Checking if the count is 0
        print("Accession number not found.")                                            # Printing a message
        return None                                                                     # Returning None

    gi_number = record["IdList"][0]                                                     # Getting the GI number
    handle = Entrez.efetch(db="protein", id=gi_number, rettype="gb", retmode="text")    # Defining the handle
    record = SeqIO.read(handle, "genbank")                                              # Reading the handle
    handle.close()                                                                      # Closing the handle

    accession_number = record.annotations["accessions"][0]                              # Getting the accession number
    print("Accession number:", accession_number)                                        # Printing the accession number
    return accession_number                                                             # Returning the accession number

# Main execution
accession_number = get_accession_number()                                               # Getting the accession number

# Download protein from NCBI 
def download_protein_sequence(accession):                                   # Defining the function
    Entrez.email = 'jungwirth@campus.tu-berlin.de'                          # Set your email address for Entrez
    handle = Entrez.efetch(db='protein',                                    # Defining the handle
                           id=accession,                                    # Defining the accession number
                           rettype='fasta',                                 # Defining the return type
                           retmode='text')                                  # Defining the return mode
    return handle.read()                                                    # Returning the record

# Download core-MetE protein sequence
accession = accession_number                                                # Defining the accession number
sequence = download_protein_sequence(accession)                             # Defining the sequence
print(sequence)                                                             # Printing the sequence

#Save sequence as a fasta file  #needs the SeqIO module
file_path = "/Users/lindajungwirth/python/GeoBT/sequence.fasta"             # Defining the file path
with open(file_path, "w") as file:                                          # Opening the file
    file.write(sequence)                                                    # Writing the sequence to the file
print("Saved sequence to", file_path)                                       # Printing the file path
