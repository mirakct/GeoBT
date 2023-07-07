from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Align #von ChatGPT & Copilot wurde noch vorgeschlagen, pairwise2 zu benutzen, aber das ist deprecated --> wird nicht mehr weiterentwickelt
from Bio import AlignIO

Entrez.email = ""

# Accession number of the protein sequence
accession_number = "WP_278812783"

# Perform BLAST search
result_handle = NCBIWWW.qblast("blastp", "nr", accession_number)

# Parse the BLAST result
blast_record = NCBIXML.read(result_handle)

print("Number of hits:", len(blast_record.alignments))

# Access the first hit (most significant)
alignment = blast_record.alignments[0]

# Download the sequence of the hit
hit_accession = alignment.accession
fasta_handle = NCBIWWW.qblast("blastp", "nr", hit_accession)
fasta_record = fasta_handle.read()

# Perform pairwise global alignment with the downloaded sequence
seq1 = blast_record.query
seq2 = fasta_record

aligner = Align.PairwiseAligner()
#alignments = aligner.align([seq1], [seq2])
alignments = aligner.align(seq1, seq2, strand="+")

# Print the alignments and alignment characteristics
for alignment in alignments:
    print(alignment)
    print("Alignment length:", alignments[0].aligned[0].end - alignments[0].aligned[0].start)
    print("Alignment Score:", alignment.score)
    print("Identity:", alignment.identities / alignment.length)
    print("Similarity:", alignment.similarity / alignment.length)
    print("Gaps:", alignment.gaps / alignment.length)
    print("---")
    print("Length:", alignment.length)



# Close the handles
result_handle.close()
fasta_handle.close()
