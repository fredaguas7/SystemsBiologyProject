import cobra
from cobra.flux_analysis.variability import find_essential_genes
from cobra.flux_analysis.variability import find_essential_reactions


model = cobra.io.load_json_model("iEK1011_m7H10_media.json")

#Checking if the model is ok
solution = model.optimize()

print(model.summary())

essential_genes = find_essential_genes(model)
# Print essential genes
for gene in essential_genes:
    print(gene)

# Identify essential reactions
essential_reactions = find_essential_reactions(model)

# Print essential reactions
for reaction in essential_reactions:
    print(reaction.id, reaction.name)


essential_gene_ids = [gene.id for gene in essential_genes]

def get_fasta_sequences(fasta_file, essential_gene_ids):
    fasta_sequences = []

    with open(fasta_file, "r") as f:
        lines = f.readlines()

    for target_id in essential_gene_ids:
        truncated_target_id = target_id[:-1]
        fasta_sequence = ""
        for i, line in enumerate(lines):
            if line.startswith(">") and truncated_target_id in line:
                for j in range(i + 1, len(lines)):
                    if lines[j].startswith(">"):
                        break
                    fasta_sequence += lines[j].strip()
                break
        if fasta_sequence:
            fasta_sequences.append((target_id, fasta_sequence))

    return fasta_sequences

fasta_file = "uniprot-compressed_true_download_true_format_fasta_query__28_28prote-2023.04.22-13.59.24.74.fasta"
# Call the function to get the matched FASTA sequences
fasta_sequences = get_fasta_sequences(fasta_file, essential_gene_ids)

# Save the matched FASTA sequences to a new file
with open("matched_sequences.fasta", "w") as output_file:
    for gene, sequence in fasta_sequences:
        output_file.write(f"> {gene}\n{sequence}\n")



fasta_file = "uniprot-compressed_true_download_true_format_fasta_query__28_28prote-2023.04.22-13.59.24.74.fasta"

fasta_sequences = get_fasta_sequences(fasta_file, essential_gene_ids)

for gene, sequence in fasta_sequences:
    print(f"> {gene}\n{sequence}\n")

################################################################################
##### Run Omega Fold and P2Rank outside python before running the rest #########
################################################################################

from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd

#Separate the fasta file in multiple fastas


# Define the input and output filenames
input_file = "matched_sequences.fasta"

ids = []

# Open the input FASTA file and loop over each record
for record in SeqIO.parse(input_file, "fasta"):
    # Create an output filename using the sequence identifier
    output_file = f"{record.id}.fasta"
    # Write the single sequence to the output file
    with open(output_file, "w") as out_handle:
        SeqIO.write(record, out_handle, "fasta")

    ids.append(record.id)

#Blast for all of them
    
for i in ids:
    
    my_query = SeqIO.read(i+".fasta", format="fasta")
    result_handle = NCBIWWW.qblast(program="blastp", database="refseq_protein", sequence = my_query.seq, entrez_query = "human[ORGN]")
    blast_result = open(i+"_blast.xml", "w")
    blast_result.write(result_handle.read())
    blast_result.close()
    result_handle.close()
    


#Create 0 hits IDs

zero_hits = []

for i in ids:
    
    # Define the input XML filename
    input_file = i + "_blast.xml"
    
    # Open the XML file and parse the results
    with open(input_file, "r") as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        # Loop over each BLAST record and count the hits
        hit_count = 0
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                hit_count += 1

    if hit_count == 0:
        zero_hits.append(i)

# Rank P2Rank Outputs

scores = pd.DataFrame({
    "Score": [],
    "Probability": []
}, index=[])

for i in zero_hits:
    path = 'C:/Users/freda/Desktop/MCBBI/BS/Project/test_output/predict_' + i[2:6]
        
    df = pd.read_csv(path + "/" + i[2:6] + '_predictions.csv')
    try:
        score = df['   score'][0]
        probability = df[' probability'][0]
    
        scores.loc[i] = [score, probability]
    except:
        scores.loc[i] = [0, 0]

scores.sort_values(by=['Score'], ascending = False)
    
scores.head()    

# Check if the knockout of identified targets really yields 0 biomass

top = ["Rv1391","Rv1631","Rv2443","Rv3215","Rv3818"]

print('complete model: ', model.optimize().objective_value)
with model:
    model.genes.Rv1391.knock_out()
    print( 'Rv1391 knocked out: ', model.optimize().objective_value)
    
with model:
    model.genes.Rv1631.knock_out()
    print( 'Rv1631 knocked out: ', model.optimize().objective_value)
    
with model:
    model.genes.Rv2443.knock_out()
    print( 'Rv2443 knocked out: ', model.optimize().objective_value)

with model:
    model.genes.Rv3215.knock_out()
    print( 'Rv3215 knocked out: ', model.optimize().objective_value)

with model:
    model.genes.Rv3818.knock_out()
    print( 'Rv3818 knocked out: ', model.optimize().objective_value)
