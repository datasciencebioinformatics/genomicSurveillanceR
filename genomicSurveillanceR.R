#########################################################################
# https://github.com/datasciencebioinformatics/genomicSurveillanceR
# genomicSurveillanceR.R
# Script to process fasta file and metadata from GISAID data base
# Inputs  --
# Input 1 : .fasta sequence file file
# Input 2 : .tsv patient metadata file
##########################################################################
# Load libraries to be used in the pipline
library("base")         
library("caret")
library("dplyr")
library("factoextra")
library("ggfortify")
library("ggnewscale")
library("ggplot2")
library("ggtree")
library("gridExtra")
library("multiplex")
library("magrittr")
library("RColorBrewer")
library("scales")
library("stringr")
library("tidytree")
library("treeio")
library("viridis")
library("viridis")                 
library("viridisLite")
library("xlsx")
library("msa")
#########################################################################
# Path to fasta file as retrieved from GISAID
fasta_path="./1684404410482.sequences.fasta"

# Path to TSV file as retrieved from GISAID
# Years 2022-2023
metadata_path="./1684404410482.metadata.tsv"

##########################################################################
# Data is submitted to nextclade as downloade:
# 1) Sequences
##########################################################################
# First load metadata file
metadata_table=read.table(file = metadata_path, sep = '\t', header = TRUE)

# Set rownames of metadata
rownames(metadata_table)<-metadata_table$gisaid_epi_isl
##########################################################################
# Check columns to be used in the analsis
# Adjust tsv metadata accordinly
head(metadata_table)[,1] # strain             ex.hCoV-19/Brazil/SP-NVBS26147GENOV852200493362/2022
head(metadata_table)[,2] # virus              ex.betacoronavirus
head(metadata_table)[,3] # gisaid_epi_isl     ex.gisaid_epi_isl	
head(metadata_table)[,5] # date               ex.2022-11-14

head(metadata_table)[,6]  # region             ex.South America
head(metadata_table)[,4] # genbank_accession  ex.?
head(metadata_table)[,7]  # country            ex.Brazil
head(metadata_table)[,8]  # division           ex.Sao Paulo	
head(metadata_table)[,9]  # location           ex.
head(metadata_table)[,10] # region_exposure    ex.South America

head(metadata_table)[,11]  # country_exposure  ex.Brazil
head(metadata_table)[,12]  # division_exposure ex.Sao Paulo
head(metadata_table)[,13]  # segment           ex.genome
head(metadata_table)[,14]  # length            ex.29765
head(metadata_table)[,15]  # host              ex.Human

head(metadata_table)[,16]  # age               ex.24
head(metadata_table)[,17]  # sex               ex.Female
head(metadata_table)[,18]  # Nextstrain_clade  ex.?
head(metadata_table)[,19]  # pangolin_lineage  ex.BE.10
head(metadata_table)[,20]  # GISAID_clade      ex.GRA

head(metadata_table)[,21]  # originating_lab  ex.DASA
head(metadata_table)[,22]  # submitting_lab   ex.DASA
head(metadata_table)[,23]  # authors          ex.Felipe
head(metadata_table)[,24]  # url              ex.https://www.gisaid.org/
head(metadata_table)[,25]  # title            ex.?

head(metadata_table)[,21]  # originating_lab  ex.DASA
head(metadata_table)[,22]  # submitting_lab   ex.DASA
head(metadata_table)[,23]  # authors          ex.Felipe
head(metadata_table)[,24]  # url              ex.https://www.gisaid.org/
head(metadata_table)[,25]  # title            ex.?

head(metadata_table)[,26]  # paper_url              ex.?
head(metadata_table)[,27]  # date_submitted         ex.2022-12-27
head(metadata_table)[,28]  # purpose_of_sequencing  ex.?
##########################################################################
# Make a copy of the table sorted by date 
metadata_table_sub<-metadata_table[order(metadata_table$date, decreasing = TRUE),]

# Obtain the day from the date
metadata_table_sub$Day<-as.numeric(substring(metadata_table$date,9,10))

# Obtain the year
metadata_table_sub$year<-as.numeric(substring(metadata_table_sub$date,1,4))

# Obtain the month from the date
metadata_table_sub$Month<-substring(metadata_table_sub$date,6,7)

# Create month-year field
metadata_table_sub$MonthYear=paste(metadata_table_sub$Month,metadata_table_sub$year,sep="-")

#########################################################################  
# Arrange entries according to the date (sort table)
metadata_table_sub=metadata_table_sub[order(as.Date(paste("01",metadata_table_sub$MonthYear,sep="-"), "%d-%m-%Y")),]

#########################################################################  
# Checkpoint : plot number of cases per age
# Count the number of cases per sex and age
metadata_table_sub_age=data.frame(table(metadata_table_sub[,c("age","sex")])[,c("Female","Male")])

# Sort the table by the number of cases
metadata_table_sub_age<-metadata_table_sub_age[order(-metadata_table_sub_age$Freq),]
	
# Remove columns with unknown age
metadata_table_sub_age<-metadata_table_sub_age[metadata_table_sub_age$age!="?",]

# If contains word unknown - set to na
metadata_table_sub_age[metadata_table_sub_age$age=="unknown","age"]=NA

# If less than 1 year - set to zero years
metadata_table_sub_age[grepl("month", metadata_table_sub_age$age, fixed=TRUE),"age"]=0

# If less than 1 year - set to zero years
metadata_table_sub_age[grepl("m", metadata_table_sub_age$age, fixed=TRUE),"age"]=0

# If less than 1 month - set to zero years
metadata_table_sub_age[grepl("day", metadata_table_sub_age$age, fixed=TRUE),"age"]=0

# Convert age to numeric
metadata_table_sub_age$age<-as.numeric(as.vector(metadata_table_sub_age$age))

# Calculate categories for numeric age - from 0 to 105 by groups 
categories=cut(as.numeric(as.vector(metadata_table_sub_age$age)), seq(0, 105, 10))

# Set levels for column categories
levels=levels(cut(as.numeric(as.vector(metadata_table_sub_age$age)), seq(0, 105, 10)))

# Adjust string label
categories_replaced=gsub("\\]", " anos", gsub("\\(", "", gsub(",", "-", categories)))

# Adjust string label
levels_replaced=gsub("\\]", " anos", gsub("\\(", "", gsub(",", "-", levels)))

# Add column categorie to metadata_table_sub_year
metadata_table_sub_age$Cat=categories_replaced

# Replace categories for age = 0
metadata_table_sub_age[which(metadata_table_sub_age$age==0),"Cat"]<-"0-10 anos"

# Remove entries with NA
metadata_table_sub_age=metadata_table_sub_age[complete.cases(metadata_table_sub_age),]

# Rename columns from frequêncy to cases
colnames(metadata_table_sub_age)[3]<-"Cases"

# Replace age data to category
metadata_table_sub_age$age<-metadata_table_sub_age$Cat

# Adjust vector for sex
metadata_table_sub_age$sex=as.vector(metadata_table_sub_age$sex)

# Create factor from vector
metadata_table_sub_age$sex=factor(metadata_table_sub_age$sex)
  
# Set a vector for the breaks
sequence_breaks=seq(-max(metadata_table_sub_age[metadata_table_sub_age$sex=="Female","Cases"]),max(metadata_table_sub_age[metadata_table_sub_age$sex=="Male","Cases"]),by=3)

# Create thple for "Número de casos por faixa etária"
p1<-ggplot(metadata_table_sub_age, aes(x = age, y = Cases, fill = sex)) +
  geom_bar(data = subset(metadata_table_sub_age, sex == "Female"), aes(y = -Cases), stat = "identity", position=position_dodge(width=2.5)) +
  geom_bar(data = subset(metadata_table_sub_age, sex == "Male"), stat = "identity", position=position_dodge(width=2.5))  +
  coord_flip()+ theme_bw()+ theme(legend.position="bottom",text = element_text(size = 12))  +ggtitle("Number of cases per age group") + ylab("Number of cases") + xlab("age group")   + scale_y_continuous(breaks=sequence_breaks,labels=abs(sequence_breaks)) 
  
# Set colours   
p1<- p1  + scale_fill_viridis_d() + guides(fill=guide_legend(title="Sex")) +  theme_bw()+ theme(legend.position="right",text = element_text(size = 8))

# Save figure
ggsave(filename="./plot_number_of_cases_per_age_group.png", plot = p1, width = 8, height = 5)
#########################################################################  
# Create count table for Cases  sex and weeks
metadata_table_sub_date=data.frame(table(metadata_table_sub[,c("MonthYear","sex")])[,c("Female","Male")])

# Rename the Freq column to Cases
colnames(metadata_table_sub_date)[3]<-"Cases" 

# Set the breaks
sequence_breaks=seq(-max(metadata_table_sub_date[metadata_table_sub_date$sex=="Male","Cases"]),max(metadata_table_sub_date[metadata_table_sub_date$sex=="Female","Cases"]),by=5)

# Re-level facor 
metadata_table_sub_date$MonthYear=factor(metadata_table_sub_date$MonthYear,levels=unique(metadata_table_sub$MonthYear))

# Create the plot for "Número de casos por périodo"
p2<-ggplot(metadata_table_sub_date, aes(x = MonthYear, y = Cases, fill = sex)) +
  geom_bar(data = subset(metadata_table_sub_date, sex == "Female"), aes(y = -Cases), stat = "identity") +
  geom_bar(data = subset(metadata_table_sub_date, sex == "Male"), stat = "identity") +
  coord_flip()+ theme_bw() + scale_fill_manual(values=c(Male = "darkblue", Female = "violet")) + theme(legend.position="right",text = element_text(size = 8))+ggtitle("Number of cases per period") + xlab("date") + ylab("Number of cases")  + theme(panel.spacing = unit(0, "lines")) + scale_fill_viridis_d()+ guides(fill=guide_legend(title="Sex"))

# Save figure
ggsave(filename="./plot_number_of_cases_per_period.png", plot = p2, width = 8, height = 5)
#########################################################################      
# Create count table for Cases  city and sex
metadata_table_sub_city=data.frame(table(metadata_table_sub[,c("location","sex")])[,c("Female","Male")])

# Sort columns by Frequency
metadata_table_sub_city<-metadata_table_sub_city[order(-metadata_table_sub_city$Freq),]

# Create count table
metadata_table_sub_city<-metadata_table_sub_city[which(metadata_table_sub_city$location %in% metadata_table_sub_city[which(metadata_table_sub_city$sex=="Male"),"location"]),]

# Rename the Freq column to Cases
colnames(metadata_table_sub_city)[3]<-"Cases"

# a hack, but it works
all_city = as.vector(unique(metadata_table_sub_city$location))                             # Take the name of all the cities

# Create a factor fot the cities in alphabetical order
metadata_table_sub_city$order = factor(metadata_table_sub_city$location, levels =all_city) 

# Divide the count per city by 100
metadata_table_sub_city$order = as.numeric(metadata_table_sub_city$location)/100           

# Select top cities
TOP_N=10
top_cities<-as.vector(metadata_table_sub_city[metadata_table_sub_city$sex=="Male","location"][1:TOP_N])

# Create a column to register the label
metadata_table_sub_city$Labels=""

# Set labels 
metadata_table_sub_city[which(metadata_table_sub_city$location %in% top_cities),"Labels"]<-as.vector(metadata_table_sub_city[which(metadata_table_sub_city$location %in% top_cities),"location"])

# Subset top 10 entries 
metadata_table_sub_city=metadata_table_sub_city[metadata_table_sub_city$Labels!="",]

# Create factors for the labels
metadata_table_sub_city$location <- factor(metadata_table_sub_city$location, levels = unique(metadata_table_sub_city$location))

# Rename the colum Freq to Cases
colnames(metadata_table_sub_city)[3]<-"Cases"

# Compute the breaks
sequence_breaks=seq(-max(metadata_table_sub_city[metadata_table_sub_city$sex=="Female","Cases"]),max(metadata_table_sub_city[metadata_table_sub_city$sex=="Male","Cases"]),by=25)

# Create the plot for "Número de casos por cidade"
p3<-ggplot(metadata_table_sub_city, aes(x = location, y = Cases, fill = sex)) +
  geom_bar(data = subset(metadata_table_sub_city, sex == "Female"), aes(y = -Cases), stat = "identity", position=position_dodge(width=2.5)) +
  geom_bar(data = subset(metadata_table_sub_city, sex == "Male"), stat = "identity", position=position_dodge(width=2.5))  +
  coord_flip()+ theme_bw()+ theme(legend.position="bottom",text = element_text(size = 9))  +ggtitle("Number of cases per city") + ylab("Number of cases") + xlab("city")   + scale_y_continuous(breaks=sequence_breaks,labels=abs(sequence_breaks))  + guides(fill=guide_legend(title="Sex"))+ theme(legend.position="right",text = element_text(size = 8))
p3<- p3  + scale_fill_viridis_d() 

# Save figure
ggsave(filename="./plot_Numero de_casos_por_cidade.png", plot = p3, width = 8.0, height = 5.0)

##########################################################################
# Plot the frequency distribution of lineages - montly collection dates (adjust to weekly or peridodically if needed)
# Re-level facor 
metadata_table_sub$MonthYear=factor(metadata_table_sub$MonthYear,levels=unique(metadata_table_sub$MonthYear))

# Frequência temporal das pangolin_lineage
pl <- ggplot(data = metadata_table_sub,aes(x= MonthYear, fill = pangolin_lineage, colour = pangolin_lineage))
pl <- pl + geom_bar(stat="count", position = "fill")
pl <- pl  + theme(axis.text.x = element_text(angle = 90,hjust =0 ))
pl <- pl + scale_y_continuous(labels = scales::percent)+ theme_bw() + theme(axis.text.x = element_text(angle = 90,hjust =0 ))  + theme(legend.position="b",text = element_text(size = 12))
pl1<-pl + theme(legend.position="right",legend.key.size = unit(0.2, 'cm'))+ggtitle("Temporal frequency of lineages") + xlab("data") + ylab("Frequência") 

# Save figure
ggsave(filename="./plot_temporal_frequency_of_lineages.png", plot = pl1, width = 10, height = 5)

#########################################################################
# Plot the distribution alines of lineages - montly collection dates (adjust to weekly or peridodically if needed)
# Calculate count tables
metadata_table_sub_lineage_counts=table(metadata_table_sub$pangolin_lineage,metadata_table_sub$MonthYear)

# Normalize count data
normalized_lineage_counts <- preProcess(data.frame(metadata_table_sub_lineage_counts), method=c("range"))
normalized_lineage_counts <- predict(normalized_lineage_counts, data.frame(metadata_table_sub_lineage_counts))

# Renames collumns
colnames(normalized_lineage_counts)<-c("lineages","date","Freq")

# stacked area chart
pl2 <- ggplot(normalized_lineage_counts, aes(x=date, y=Freq, group=lineages, fill=lineages)) + geom_area() + theme(axis.text.x = element_text(angle = 90,hjust =0 ))  + theme(legend.position="b",text = element_text(size = 12)) + theme_bw() + theme_bw() + theme(axis.text.x = element_text(angle = 90,hjust =0 ))  
    
# stacked area chart
pl2<-pl2 + theme(legend.position="right",legend.key.size = unit(0.2, 'cm'))+ggtitle("Scheme of lineage distribution") + xlab("data") 

# Save figure
ggsave(filename="./plot_scheme_of_lineage_distribution.png", plot = pl2, width = 10, height = 5)

#########################################################################
# Calculate the phylogeny from tree and alignment
# Load tree
myTree <- ape::read.tree(tree_path)

# Select data set
df=data.frame(lineage=metadata_table_sub$pangolin_lineage,id=metadata_table_sub$strain)

# Set rownames
rownames(df)<-metadata_table_sub$strain

# Set colnames
colnames(df)<-c("lineages","id")

# Set colnames
df=data.frame(df[myTree$tip.label,])

# Create ggtree plots
circ <- ggtree(myTree,branch.length='none', layout = "circular") + theme_tree2() 

# gheatmap
p1 <- gheatmap(circ, df,offset = 2, width=0.25) + scale_fill_viridis_d() + ggtitle("Lineages phylogeny")

# Save plots
ggsave(filename="./Phylogeny1.png", plot = p1, width=5, height=5)  
#######################################################################
library(ggfortify)
# Plot the principal component analysis
# Calculate the phylogeny from tree and alignment
# Calculate the clusters
cluster_myTree<-hclust(d = dist(cophenetic(myTree)))

# Calculate distance tree
d = dist(cophenetic(myTree))
pca_res <- prcomp(d, scale. = TRUE)
pca_plot=autoplot(pca_res, data =df, colour="lineages")
pca_plot = pca_plot +  theme_bw() +ggtitle("Lineages PCA")+ scale_fill_viridis_d() 
ggsave(filename="./PCA_lineages.png", plot = pca_plot, width=8, height=4) 

#######################################################################
# Calculate outliers
# First compute PCA
pca <- prcomp(d, scale. = TRUE, rank. = 10)
U <- pca$x

U2 <- U
U2[1, 1] <- 30
U3=U2

# Compute outliers
outilers=apply(U3, 2, function(x) which( abs(x - mean(x)) > (6 * sd(x)) )) %>%
  Reduce(union, .)
U3=data.frame(U3)
U3$Label=""
U3$Outliers=""  
U3$colour="black"	  
U3[outilers,"colour"]="red"
U3[outilers,"Outliers"]=rownames(U3)[outilers]
 
# Plot ouuliers 
Plot_outliers=qplot(PC1, PC2, data = U3, colour = colour, label=U3$Outliers) + theme_bw() + scale_color_manual(values=c("black", "red")) + geom_text(vjust = 0, nudge_y = 0.5) + ggtitle("Possible new variants")

# ggsave
ggsave(filename="./PCA_outliners.png", plot = Plot_outliers, width=8, height=4) 
