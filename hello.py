#download a sequence from ncbi and save it as a fasta file
Entrez.email = "linda.jungwirth1@gmail.com"
handle = Entrez.efetch(db="nucleotide", id="NM_002299.3", rettype="fasta", retmode="text")
seq_record = SeqIO.read(handle, "fasta")
handle.close()
print(seq_record.id)
print(repr(seq_record.seq))
print(len(seq_record))
SeqIO.write(seq_record, "NM_002299.fasta", "fasta")
