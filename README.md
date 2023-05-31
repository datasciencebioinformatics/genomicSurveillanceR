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

Outliners analises (possible new variant)

### Prerequisites
sudo apt-get install mafft

sudo apt-get install fasttree

sudo apt-get install r-base

R libraries are described in the .r file

## Pre-processing
### Input and output
fasta_path="./1684404410482.sequences.fasta"

aln_path="./1684404410482.sequences.aln"

tree_path="./1684404410482.sequences.tree"

### Alignment
mafft --auto $fasta_path > $aln_path

### Tree
fasttree -nt $aln_path > $tree_path

## Usage
First, execute the pre-processing steps to compute alignement (mafft) and calculate tree (fasttree). Second, follow genomicSurveillanceR.script steps.

