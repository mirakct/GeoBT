from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Align
from Bio import AlignIO
import datetime

import tkinter as tk #hello

Entrez.email = "karim.salama@campus.tu-berlin.de"


#a function to get the protein sequence from the accession number
def get_protein_sequence(accession_number):
    handle = Entrez.efetch(db="protein", id=accession_number, rettype="fasta", retmode="text")
    record = handle.read()
    handle.close()
    return record

#a function to blast the protein accession number against the nr database


def blast_protein_sequence(protein_sequence, hl=5, entrez='none', threshold='none'):
    try:
        result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence,hitlist_size=hl,filter=entrez,i_thresh=threshold)
        with open('blast.xml', 'w') as save_file:
            save_file.write(result_handle.read())
        result_handle.close()
        print("Blast done")
    except:
        errorout = open('errorlog.txt','a')
        errorout.write( str(datetime.datetime.now) +':' + 'failed blast database query' + '\n')
        errorout.close()
        return

    out=[[0 for _ in range(13)] for _ in range(hl)]

    with open('blast.xml', 'r') as blast_file:
        blast_record = NCBIXML.parse(blast_file)        

        i=0
        for alignments in blast_record:
            
            for alignment in alignments.alignments:
                gi,acc,_ = alignment.hit_id.split('|')
                prt,org = alignment.hit_def.split('[',1)
                org,_= org.split(']',1)
                out[i][0] =  acc                      #accession number
                out[i][1] =  prt                      #protein name
                out[i][2] =  org[:-1]                 #organism
                out[i][3] =  alignment.length         #length of the protein

                for hsp in alignment.hsps:
                    out[i][4] = hsp.expect               #e-value
                    out[i][5] = hsp.align_length         #alignment length
                    out[i][6] = hsp.score                #score
                    out[i][7] = hsp.identities           #identities
                    out[i][8] = hsp.positives            #positives
                    out[i][9] = hsp.gaps                 #gaps

                    out[i][10] =  hsp.query              #query
                    out[i][11] =  hsp.match              #match
                    out[i][12] =  hsp.sbjct              #subject
                i+=1
    return out


# on install change changge line_width in __init__.py of Bio.Align to 99999999 #####################
# alignment using PairwiseAligner using substitution matrix BLOSUM62
def align_sequences(seq1, seq2,scope='global'):
    
    aligner = Align.PairwiseAligner( mode=scope)
    aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    alignments = aligner.align(seq1, seq2)
    alignments = list(alignments)
    alignment = alignments[0]
    print("Alignment done")
    scores = alignment.score
    #len = len(alignments)
    print("Score: ", scores)
    print(alignment)
    return alignment,scores





