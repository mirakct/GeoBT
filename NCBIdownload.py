# Import important packages
from Bio import Entrez
from Bio import SeqIO   

# Provide your email address to NCBI
Entrez.email = "linda.jungwirth1@gmail.com"

#Get accession number
def get_accession_number():
    search_term = "Dehalococcoides mccartyi[Organism] AND methionine synthase[Protein]"
    handle = Entrez.esearch(db="protein", term=search_term, retmax=1)
    record = Entrez.read(handle)
    handle.close()

    if int(record["Count"]) == 0:
        print("Accession number not found.")
        return None

    gi_number = record["IdList"][0]
    handle = Entrez.efetch(db="protein", id=gi_number, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    accession_number = record.annotations["accessions"][0]
    print("Accession number:", accession_number)
    return accession_number

# Main execution
accession_number = get_accession_number()

# Download protein from NCBI 
def download_protein_sequence(accession):                                   # Defining the function
    Entrez.email = 'linda.jungwirth1@gmail.com'                             # Set your email address for Entrez
    handle = Entrez.efetch(db='protein',                                    # Defining the handle
                           id=accession,                                    # Defining the accession number
                           rettype='fasta',                                 # Defining the return type
                           retmode='text')                                  # Defining the return mode
    return handle.read()                                                    # Returning the record

# Download core-MetE protein sequence
accession = accession_number                                                  # Defining the accession number
sequence = download_protein_sequence(accession)                             # Defining the sequence
print(sequence)                                                             # Printing the sequence

#Save sequence as a fasta file  #needs the SeqIO module
file_path = "/Users/lindajungwirth/python/sequence.fasta"                   # Defining the file path
with open(file_path, "w") as file:                                          # Opening the file
    file.write(sequence)                                                    # Writing the sequence to the file
print("Saved sequence to", file_path)                                       # Printing the file path