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

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEhBgKhGgBUhVs9lUASPKc7Y7FocSmb9Pc2wU5n7Cn4Gfu9bxC5w-DVGpvFpzPYSVVHpnqpSm5AYAgY_d1MX45jGJxbMp5c2dAUbAdyj7xfTPBI0Lp0Dr-5ce-mnrXqaseWkd3ISv7JavAtfzqBvKmP3Z3bYgPGBUbsGXI0wna_cqXsPN-64lvZJPRkh/s2400/plot_number_of_cases_per_period.png)

### Number of cases per age group and per sex (male/female)

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEhgv7y0jBPg2idRTlitSnGjzQf7WO-xOg66QlnvOa8wBRi0PuTyekqSpjQEHr7ze4ra03N3g2t-6N7X2Q8MGJfjtgEepL6NTlUh-UhFkmLMY_XQz634SryUaEw7TkTaEiutdZb1NWv1PdJ5qYFrZ1vJMOdgRmtN7D0J-5w5Xnt0zeFeiFaTNXvqQfFR/s2400/plot_number_of_cases_per_age_group.png)

### Number of cases per city and per sex (male/female)

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEjK2oxEZS-3TTyXxXCkobNhZD7FkA9CyBCZrD1vuZquiSw-Mak3BjJaTir2rjR4zbV_jBrNpxC_KJR5CwT8akWkz1YlPaMVAt6XoJxADPWJMbIkgUBXWUdpmM9LsvEqaqhnkzulq7-04EBNKtbhPtcwVjSPfykHx9ZfKtlDlYLWIQTWoUuQbFkl4zPy/s2400/plot_Numero%20de_casos_por_cidade.png)

### Temporal distribution of lineages

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEja2vS2rgCdAnwyOgzDPt1vPkGv2cwVJVXPHRFhInu5GNx7he02_VVQ6jcgnu8ZgOaxbZo_lPQut3-JKq46JZbobGasSyuBqnjPZsqdMt5EQ8dEd_RcD7e0gQtx98RByNXbAtTpNvJ6EN0La400vOpSfeFnMgU-uE1ERSwz8P0Xwvv_798PA0-atrpL/s3000/plot_temporal_frequency_of_lineages.png)

### Scheme for frequency of lineages

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEjhSFf6-usViD9su19Fmpu0h1kj-f-ZuZ_JCpodF3DaMFea4VQr_4Nen9TuUY57TEARLYg8sIS6amWQzvCyO0rxwJ5Zft6OYmpf1oN7_qJtLXNvn4A7oTjzStW9nHrmoHlGN8Fn-XjXdwiws4JS0jkg1ayokqpu5TVsvATl5hlQciBybu6a7zGNCkWQ/s3000/plot_scheme_of_lineage_distribution.png)

### Principal component of samples per lineage

![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEiwOwkFtnVUK5QNCzPX_kOlf_p50Xv1sMwOTdOykSHcu0W6Fd8_iyBVk5qVGNPKH5bVGwXMAZZwGoFtHyvB_gYGDIpoBPCWcwe-XGY5b4plW7KJJXMMyYCJngZq1RDLKr5wOfbtKD7RDr7Qv0ZtRu6ZoiPQs4TkeI-kwbY-hJ39K_FeJGx3HI3yWsKg/s2400/PCA_lineages.png)

### Outlier analyses (possible new variant)![](https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEhbR0ffb8JRs3ArB3armbBuTqZOZEABcIj9Bcl603JshKJ5nd6blVcVvR6dZ7hmXv2MayoO0SBCR0cakyZmQVKeKT5vfV_AaZhMa5ErlN0feVxU6lI7mllqZXRpa1vnNc8sLTWxlAzYDnPBsioSdM7XJw9zmR4bisxn4kkYq3Tm0ZXyWzpPCiBDiDgN/s2400/PCA_outliners.png))
