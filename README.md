<!-- GETTING STARTED -->
# genomicSurveillanceR 
Script for bioinformatics of genomic surveillance from GISAID database. In the example, GISAID data from the Brazilian state of Sao Paulo (11-2022 to 04-2023) will be used.

###  Analysis
Number of cases per week and per sex (male/female) 

Number of cases per age group and per sex (male/female)

Number of cases per city and per sex (male/female)

Temporal distribution of lineages

Temporal frequence of lineages

Principal component of samples per lineage

Outliners analises (possible new variants)

## Pre-processing
### Set input (fasta) and output (.aln and .tree ) files
fasta_path="./1684404410482.sequences.fasta"

aln_path="./1684404410482.sequences.aln"

tree_path="./1684404410482.sequences.tree"

### Generate alignment
mafft --auto $fasta_path > $aln_path

### Calculate tree
fasttree -nt $aln_path > $tree_path

## Usage
First, execute the pre-processing steps to compute alignement (mafft) and calculate tree (fasttree).

