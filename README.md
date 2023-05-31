# genomicSurveillanceR
Script to process fasta file and metadata from GISAID data base

# Analysis
Number of cases per week and per sex (male/female)
 
Number of cases per age group and per sex (male/female)

Number of cases per city and per sex (male/female)

Temporal distribution of lineages

Temporal frequence of lineages

Principal component of samples per lineage

Outliners analises (possible new variants)

# Pre-processing
# Path to fasta file as retrieved from GISAID
fasta_path="./1684404410482.sequences.fasta"

# Obs. Use mafft and fasttree to calculate alignment and generate tree
# Aln path
aln_path="./1684404410482.sequences.aln"

# Aln path
tree_path="./1684404410482.sequences.aln"

# Run nextstrain online to obtain 
# Aligned by nextrain
# Generate alignment
mafft --auto $fasta_path > $aln_path

# Here I generate the tree with iqtree
fasttree -nt $aln_path > $tree_path



