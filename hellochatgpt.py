from Bio import SeqIO
from Bio import Entrez

Entrez.email = 'your_email@example.com'  # Replace with your email address

# Fetch the sequence from NCBI
handle = Entrez.efetch(db='nucleotide', id='NC_000913.3', rettype='fasta')
record = SeqIO.read(handle, 'fasta')
handle.close()

output_filename = 'sequence.fasta'  # Specify the output filename

# Save the sequence as a FASTA file
SeqIO.write(record, output_filename, 'fasta')
