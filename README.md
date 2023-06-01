<!-- GETTING STARTED -->

# genomicSurveillanceR 
Script for bioinformatics of genomic surveillance from GISAID database. In the example, GISAID data from the Brazilian state of Sao Paulo (11-2022 to 04-2023) will be used.

### Installation and pre-requisites
sudo apt-get install mafft

sudo apt-get install fasttree

sudo apt-get install r-base

R libraries are described in the .r file

## Pre-processing
fasta_path="./1684404410482.sequences.fasta"

metadata_path="./1684404410482.metadata.tsv"

aln_path="./1684404410482.sequences.aln"

tree_path="./1684404410482.sequences.tree"

### Alignment
mafft --auto $fasta_path > $aln_path

### Tree
fasttree -nt $aln_path > $tree_path

## Usage
First, execute the pre-processing steps to compute alignement (mafft) and calculate tree (fasttree). Second, follow genomicSurveillanceR.script steps.

#  Analysis
### Number of cases per period and per sex (male/female)

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEhAxgximn42GdZB_1iskWfZnhg8jivOvJVrwB1X4JiY2wTB3jHBfjlV-zNtPp4n97b_ZkOyRnPwbTJfR9I-gl-bo6Bq-H9nKYGhWQOs8d83FQrGsjWBmgxuzJQOMZhXiAUiUgPslY-YMrGEBNkrhlD1bCKZMJrxehOmz-b1WI4Mh89t5eY5_EXO5ydA/s2400/plot_number_of_cases_per_period.png)

### Number of cases per age group and per sex (male/female)

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEit5gzWYp-4fgKhZHoQUCdREPjI0_bfVxsmRW7NgnvGgNVmkFoTRudttv9hqccY74XfIyMP4ierENLiPVaCO0FuOSGItSKL680IjIwrxTOIcThgESYj1-X47iPzxkms-UTjJ--DzhI6n0ECn1f6Iw5J391OaD0gZXJjaD0H6z94nb7R6VSPJ8MsciOh/s2400/plot_number_of_cases_per_age_group.png)

### Number of cases per city and per sex (male/female)

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEhAxgximn42GdZB_1iskWfZnhg8jivOvJVrwB1X4JiY2wTB3jHBfjlV-zNtPp4n97b_ZkOyRnPwbTJfR9I-gl-bo6Bq-H9nKYGhWQOs8d83FQrGsjWBmgxuzJQOMZhXiAUiUgPslY-YMrGEBNkrhlD1bCKZMJrxehOmz-b1WI4Mh89t5eY5_EXO5ydA/s2400/plot_number_of_cases_per_period.png)

### Temporal distribution of lineages

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEjhSFf6-usViD9su19Fmpu0h1kj-f-ZuZ_JCpodF3DaMFea4VQr_4Nen9TuUY57TEARLYg8sIS6amWQzvCyO0rxwJ5Zft6OYmpf1oN7_qJtLXNvn4A7oTjzStW9nHrmoHlGN8Fn-XjXdwiws4JS0jkg1ayokqpu5TVsvATl5hlQciBybu6a7zGNCkWQ/s3000/plot_scheme_of_lineage_distribution.png)

### Scheme for frequency of lineages

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEjhSFf6-usViD9su19Fmpu0h1kj-f-ZuZ_JCpodF3DaMFea4VQr_4Nen9TuUY57TEARLYg8sIS6amWQzvCyO0rxwJ5Zft6OYmpf1oN7_qJtLXNvn4A7oTjzStW9nHrmoHlGN8Fn-XjXdwiws4JS0jkg1ayokqpu5TVsvATl5hlQciBybu6a7zGNCkWQ/s3000/plot_scheme_of_lineage_distribution.png)

### Principal component of samples per lineage

### Outlier analyses (possible new variant)


