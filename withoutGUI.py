from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Align, pairwise2

# Accession number of the protein sequence
accession_number = "WP_278812783"

# Perform BLAST search
result_handle = NCBIWWW.qblast("blastp", "nr", accession_number)

# Parse the BLAST result
blast_record = NCBIXML.read(result_handle)

# Access the first hit (most significant)
alignment = blast_record.alignments[0]

# Download the sequence of the hit
hit_accession = alignment.accession
fasta_handle = NCBIWWW.qblast("fasta", "nr", hit_accession)
fasta_record = fasta_handle.read()

# Perform pairwise global alignment with the downloaded sequence
seq1 = blast_record.query
seq2 = fasta_record

aligner = Align.PairwiseAligner()
alignments = aligner.align([seq1], [seq2])

# Print the alignments and alignment characteristics
for alignment in alignments:
    print(alignment)
    print("Alignment Score:", alignment.score)
    print("Length:", alignment.length)
    print("Identity:", alignment.identities / alignment.length)
    print("Similarity:", alignment.similarity / alignment.length)
    print("Gaps:", alignment.gaps / alignment.length)
    print("---")

# Close the handles
result_handle.close()
fasta_handle.close()
