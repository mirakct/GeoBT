# Import important packages 
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import matplotlib.pyplot as plt 
import matplotlib.pyplot as plt


#fasta_file = "sequence.fasta"                                               # Load the protein sequence from the FASTA file

def analyze(sequence):                                                       # Defining the function	
    #1. Amino acid composition
    #a) in percent
 
    try:                                                                        # "Try to run the following code"
        protein = ProteinAnalysis(sequence)                                     # Create a ProteinAnalysis object
        composition = protein.get_amino_acids_percent()                         # Calculate amino acid composition
        print(composition)                                                      # Print the result
    except Exception as e:                                                      # "If an exception occurs, run the following code"
        print("An error occurred:", str(e))                                     # Print the error message
    #b) in absolute numbers: 
    protein = ProteinAnalysis(sequence)                                         # Create a ProteinAnalysis object
    composition = protein.get_amino_acids_percent()                             # Calculate amino acid composition
    count = protein.count_amino_acids()                                         # Calculate amino acid count
    print(protein.count_amino_acids())                                          # Print the result

    #Visualisation
                                            
    labels = list(composition.keys())                                           # Convert composition dictionary to lists
    values = list(composition.values())                                         # Convert composition dictionary to lists
    plt.bar(labels, values)                                                     # Create a bar plot
    plt.xlabel('Amino Acid')                                                    # Add labels
    plt.ylabel('Percentage')                                                  
    plt.title('Amino Acid Composition')
    plt.show()                                                                  # Show the plot

    #2. Hydrophobicity 
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

    #3. molecular weight
    analyzer = ProteinAnalysis(sequence)                                        # Create a ProteinAnalysis object
    molecular_weight = analyzer.molecular_weight()                              # Calculate the molecular weight
    mw=f"Molecular Weight: {molecular_weight:.2f} g/mol"                    # Print the molecular weight

    #4. isoelectric point
    analyzer = ProteinAnalysis(sequence)                                        # Create a ProteinAnalysis object
    pI = analyzer.isoelectric_point()                                           # Calculate the isoelectric point
    ip=f"Isoelectric Point (pI): {pI:.2f}"                                  # Print the pI

    #5. solubility (=pI), instability
    analyzer = ProteinAnalysis(sequence)                                        # Create a ProteinAnalysis object
    solubility = analyzer.isoelectric_point()                                   # Calculate solubility = pI 
    instability_index = analyzer.instability_index()                            # Calculate instability index
    sol=f"Solubility: {solubility:.2f}"                                      # Print the results
    ins=f"Instability Index: {instability_index:.2f}"                      

    out = '\n'.join((mw,ip,sol,ins))                                                # Join the results  
    return out                                                     # Return the results

if __name__ == '__main__':                                                      # If the script is executed directly
    pass