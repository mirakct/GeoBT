# Import important packages 
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO

fasta_file = "sequence.fasta"                                               # Load the protein sequence from the FASTA file
record = SeqIO.read(fasta_file, "fasta")                                    # Read the FASTA file
sequence = str(record.seq)                                                  # Extract the sequence from the record

#1. Amino acid composition, 2 options
#b) try-except block (to handle any exceptions)
try:                                                                        # "Try to run the following code"
    protein = ProteinAnalysis(sequence)                                     # Create a ProteinAnalysis object
    composition = protein.get_amino_acids_percent()                         # Calculate amino acid composition
    print(composition)                                                      # Print the results   
except Exception as e:                                                      # "If an exception occurs, run the following code"
    print("An error occurred:", str(e))                                     # Print the error message
#a) easier: 
protein = ProteinAnalysis(sequence)                                         # Create a ProteinAnalysis object
composition = protein.get_amino_acids_percent()                             # Calculate amino acid composition
count = protein.count_amino_acids()                                         # Calculate amino acid count
print(protein.count_amino_acids())

#a) Visualisation
import matplotlib.pyplot as plt                                             
labels = list(composition.keys())                                           # Convert composition dictionary to lists
values = list(composition.values())                                         # Convert composition dictionary to lists
plt.bar(labels, values)                                                     # Create a bar plot
plt.xlabel('Amino Acid')                                                    # Add labels
plt.ylabel('Percentage')                                                    
plt.title('Amino Acid Composition')
plt.show()                                                                  # Show the plot

#2. Hydrophobicity: 2 options (which lead to different graphs) -> a) calculates the hydrophobicity directly for each residue using the predefined KD scale and plots the values against the residue numbers. b) calculates hydrophobicity scores using a sliding window approach with the KD scale, resulting in a smoothed representation of hydrophobicity scores across the sequence
##a) from the internet https://www.biob.in/2014/05/hydrophobicity-plot-using-biopython.html#:~:text=Hydrophobicity%20is%20the%20property%20of,exposed%20loops%20or%20buried%20residues.)
#from pylab import *
#for record in SeqIO.parse("sequence.fasta", "fasta"):                       # Read the FASTA file
#    id = record.id                                                          # Extract the ID
#    seq = record.seq                                                        # Extract the sequence
#    num_residues = len(record.seq)                                          # Calculate the number of residues
#kd = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,                        # Kyte-Doolittle hydrophobicity scale dictionary
#       'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
#       'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
#       'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }
#values = []                                                                 # Create an empty list
#for residue in seq:                                                         # Iterate over the sequence
# values.append(kd[residue])                                                 # Append the hydrophobicity value of the residue to the list    
#x_data = range(1, num_residues+1)                                           # Create a list of residue numbers
#plot(x_data, values, linewidth=1.0)                                         # Plot the hydrophobicity values
#axis(xmin = 1, xmax = num_residues)                                         # Set the x-axis limits      
#xlabel("Residue Number")                                                    # Add labels
#ylabel("Hydrophobicity")
##title("K&D Hydrophobicity for " + id)
#show()                                                                      # Show the plot

##b) with the help of ChatGPT (after several tries, first it was trying it without kd, then with kd method assuming it is on biopython, but then manually with kd (defining the kyte-doolittle dictionary manually))
import matplotlib.pyplot as plt
kyte_doolittle = {                                                          # Kyte-Doolittle hydrophobicity scale dictionary
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
    'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
}
analyzer = ProteinAnalysis(sequence)                                        # Create a ProteinAnalysis object
hydrophobicity_scores = analyzer.protein_scale(window=9,                    # Calculate hydrophobicity scores for each amino acid (The window parameter is set to 9, indicating that a sliding window of size 9 is used to calculate the hydrophobicity scores. This means that for each amino acid in the protein sequence, the hydrophobicity score is calculated based on the hydrophobicity values of the 9 neighboring amino acids (4 on each side, and the central amino acid itself)
                                               edge=1.0,                    # The edge parameter is set to 1.0, indicating that the edges of the protein sequence are treated as helices
                                               param_dict=kyte_doolittle)   # The param_dict parameter is set to the kyte_doolittle dictionary, indicating that the Kyte-Doolittle hydrophobicity scale is used to calculate the hydrophobicity scores
plt.figure(figsize=(10, 4))                                                 # Plot the hydrophobicity scores
plt.plot(hydrophobicity_scores)                                            
plt.xlabel("Residue Index")
plt.ylabel("Hydrophobicity Score")
plt.title("Hydrophobicity Plot")
plt.tight_layout()                                                          
plt.show()                                                                  # Show the plot

#3. GRAVY (using ChatGPT): GRAVY is the grand average of hydropathy. It is calculated by summing the hydropathy values of all the amino acids in the sequence and dividing by the number of residues in the sequence.
# GRAVY = (Σhydropathy values of all the amino acids in the sequence) / (number of residues in the sequence)
# GRAVY values range from -4.5 to +4.5. The more positive the GRAVY value, the more hydrophobic the protein is. The more negative the GRAVY value, the more hydrophilic the protein is.
# GRAVY values can be used to compare the hydrophobicity of different proteins. The protein with the higher GRAVY value is more hydrophobic than the protein with the lower GRAVY value.
analyzer = ProteinAnalysis(sequence)                                        # Create a ProteinAnalysis object
gravy_score = analyzer.gravy()                                              # Calculate the GRAVY score
print(f"GRAVY Score: {gravy_score:.2f}")                                    # Print the GRAVY score

#3. molecular weight (ChatGPT) 
analyzer = ProteinAnalysis(sequence)                                        # Create a ProteinAnalysis object
molecular_weight = analyzer.molecular_weight()                              # Calculate the molecular weight
print(f"Molecular Weight: {molecular_weight:.2f} g/mol")                    # Print the molecular weight

#4. isoelectric point (ChatGPT)
analyzer = ProteinAnalysis(sequence)                                        # Create a ProteinAnalysis object
pI = analyzer.isoelectric_point()                                           # Calculate the isoelectric point
print(f"Isoelectric Point (pI): {pI:.2f}")                                  # Print the pI

#5. solubility (=pI), instability
analyzer = ProteinAnalysis(sequence)                                        # Create a ProteinAnalysis object
solubility = analyzer.isoelectric_point()                                   # Calculate solubility = pI 
instability_index = analyzer.instability_index()                            # Calculate instability index
print(f"Solubility: {solubility:.2f}")                                      # Print the results
print(f"Instability Index: {instability_index:.2f}")