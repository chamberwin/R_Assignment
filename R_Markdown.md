---
title: "R_Assignment"
author: "H. Chamberlain"
date: '2022-03-10'
output: html_document
---
## Install/load packages
```{r}
library(tidyverse)  #For data maniupluation
library(dplyr)
install.packages("janitor")
install.packages("data.table")
library(janitor)
library(data.table)
```

## Upload raw data
```{r}
geno <- read.delim("https://raw.githubusercontent.com/chamberwin/BCB546-Spring2022/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt", header = TRUE, sep = "\t") #Download raw files from class repository for genotypes and SNPs
snp_pos <- read.delim("https://raw.githubusercontent.com/chamberwin/BCB546-Spring2022/main/assignments/UNIX_Assignment/snp_position.txt", header = TRUE, sep = "\t") 
```

## Inspect data
```{r}
ncol(geno) #Number of columns for genotype data is 986
nrow(geno) #Number of rows is 2781
dim(geno) #Dimensions of rows (2781) and columns (986)
length(geno) #Number of columns is 986
object.size(geno) #Size of file is 22681376 bytes
summary(geno) #Because of large number of columns, not super information
head(geno) #A little better than summary at seeing structure
ncol(snp_pos) #Number of columns in SNP position data is 15
nrow(snp_pos) #Number of rows in SNP position data is 983
dim(snp_pos) #Dimensions of rows (983) and columns (15)
length(snp_pos) #Number of columns is 15
object.size(snp_pos)
summary(snp_pos) #Size of file is 327392 bytes
head(snp_pos) #Look at data structure

```

## Prepare SNP position file for join
```{r}
snp_cut <- snp_pos %>% #Create prepared file for join
  select(SNP_ID, Chromosome, Position) %>% #Select only these columns
  arrange(SNP_ID) #Arrange by SNP_ID
view(snp_cut) #Check work
```

## Maize Join
```{r}
geno_maize <- geno %>%
  filter(Group %in% c("ZMMIL", "ZMMLR", "ZMMMR", "Sample_ID")) %>% # filter by maize groups only, including sample ID
  t() %>% # transpose - note it's added automatic column names (V1, V2, etc. in new header)
  as.data.frame() %>% # make sure it's a dataframe
  setDT(keep.rownames = TRUE) %>% #gives column 1 a placeholder name as well
  row_to_names(row_number = 1) %>% # sets first row as column names
  rename(SNP_ID = Sample_ID) # rename first column as SNP_ID so we can join on like column names
view(geno_maize) #double check joined by SNP_ID
join_maize <- snp_cut %>%
  inner_join(geno_maize, by = "SNP_ID") # joins both files - inner_join makes sure we take only SNP_IDs in both files
view(join_maize) #double check it joined correctly
```

## Teosinte Join
```{r}
geno_teosinte <- geno %>%
  filter(Group %in% c("ZMPBA", "ZMPIL", "ZMPJA")) %>% # filter by maize groups only, including sample ID
  t() %>% # transpose - note it's added automatic column names (V1, V2, etc. in new header)
  as.data.frame() %>% # make sure it's a dataframe
  setDT(keep.rownames = TRUE) %>% #gives column 1 a placeholder name as well
  row_to_names(row_number = 1) %>% # sets first row as column names
  rename(SNP_ID = Sample_ID) # rename first column as SNP_ID so we can join on like column names
view(geno_teosinte) #double check it is orded by SNP_ID
join_teosinte <- snp_cut %>%
  inner_join(geno_teosinte, by = "SNP_ID") # joins both files - inner_join makes sure we take only SNP_IDs in both files
view(join_teosinte) #double check it joined correctly
```

## Sort Joined Maize
```{r}
sort_maize <- join_maize %>%
  arrange(as.numeric(Chromosome), as.numeric(Position)) # sorts by chromosome number and position, in ascending order
view(sort_maize) #Check work
chr1_maize <- sort_maize %>% #Creates object for chr1
  filter(Chromosome == 1) #filters for chromosome 1 only
view(chr1_maize) #Check work
chr2_maize <- sort_maize %>% #Repeat for all chromosomes
  filter(Chromosome == 2)
view(chr2_maize)
chr3_maize <- sort_maize %>%
  filter(Chromosome == 3)
view(chr3_maize)
chr4_maize <- sort_maize %>%
  filter(Chromosome == 4)
view(chr4_maize)
chr5_maize <- sort_maize %>%
  filter(Chromosome == 5)
view(chr5_maize)
chr6_maize <- sort_maize %>%
  filter(Chromosome == 6)
view(chr6_maize)
chr7_maize <- sort_maize %>%
  filter(Chromosome == 7)
view(chr7_maize)
chr8_maize <- sort_maize %>%
  filter(Chromosome == 8)
view(chr8_maize)
chr9_maize <- sort_maize %>%
  filter(Chromosome == 9)
view(chr9_maize)
chr10_maize <- sort_maize %>%
  filter(Chromosome == 10)
view(chr10_maize)
write.table(chr1_maize, file = "chr1_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE) #creates .txt file for maize chromosome 1
write.table(chr2_maize, file = "chr2_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE) #Repeat for all chromosomes
write.table(chr3_maize, file = "chr3_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr4_maize, file = "chr4_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr5_maize, file = "chr5_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr6_maize, file = "chr6_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr7_maize, file = "chr7_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr8_maize, file = "chr8_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr9_maize, file = "chr9_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr10_maize, file = "chr10_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
```

## Sort Joined Teosinte
```{r}
sort_teosinte <- join_teosinte %>%
  arrange(as.numeric(Chromosome), as.numeric(Position)) #arrange teosinte object by chromosome and position
view(sort_teosinte) #Check work
chr1_teosinte <- sort_teosinte %>%  #create sorted teosinte objects
  filter(Chromosome == 1) #Filter by chromosome 1
view(chr1_teosinte) #Check work
chr2_teosinte <- sort_teosinte %>% #Repeat for all chromosomes
  filter(Chromosome == 2)
view(chr2_teosinte)
chr3_teosinte <- sort_teosinte %>%
  filter(Chromosome == 3)
view(chr3_teosinte)
chr4_teosinte <- sort_teosinte %>%
  filter(Chromosome == 4)
view(chr4_teosinte)
chr5_teosinte <- sort_teosinte %>%
  filter(Chromosome == 5)
view(chr5_teosinte)
chr6_teosinte <- sort_teosinte %>%
  filter(Chromosome == 6)
view(chr6_teosinte)
chr7_teosinte <- sort_teosinte %>%
  filter(Chromosome == 7)
view(chr7_teosinte)
chr8_teosinte <- sort_teosinte %>%
  filter(Chromosome == 8)
view(chr8_teosinte)
chr9_teosinte <- sort_teosinte %>%
  filter(Chromosome == 9)
view(chr9_teosinte)
chr10_teosinte <- sort_teosinte %>%
  filter(Chromosome == 10)
view(chr10_teosinte)
write.table(chr1_teosinte, file = "chr1_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE) #creates .txt files for teosinte
write.table(chr2_teosinte, file = "chr2_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr3_teosinte, file = "chr3_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr4_teosinte, file = "chr4_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr5_teosinte, file = "chr5_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr6_teosinte, file = "chr6_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr7_teosinte, file = "chr7_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr8_teosinte, file = "chr8_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr9_teosinte, file = "chr9_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr10_teosinte, file = "chr10_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
```

## Maize in descending order with ? replaced with dash
```{r}
sort_desc_maize <- join_maize %>% #Create object for sorted descending maize
  arrange(as.numeric(Chromosome), desc(as.numeric(Position))) #sort maize in descending position
view(sort_desc_maize) #Check work
maize_dash <- data.frame(lapply(sort_desc_maize, gsub, pattern = "[?]", replacement = "-")) #replace ? within brackets because it is wildcard
view(maize_dash) #Check work
chr1_maize_dash <- maize_dash %>% #creates descending sorted maize object where ? is replaced with -
  filter(Chromosome == 1) #Filter by chromosome 1
view(chr1_maize_dash) #Check work
chr2_maize_dash <- maize_dash %>% #Repeat for all chromosomes
  filter(Chromosome == 2)
view(chr2_maize_dash)
chr3_maize_dash <- maize_dash %>%
  filter(Chromosome == 3)
view(chr3_maize_dash)
chr4_maize_dash <- maize_dash %>%
  filter(Chromosome == 4)
view(chr4_maize_dash)
chr5_maize_dash <- maize_dash %>%
  filter(Chromosome == 5)
view(chr5_maize_dash)
chr6_maize_dash <- maize_dash %>%
  filter(Chromosome == 6)
view(chr6_maize_dash)
chr7_maize_dash <- maize_dash %>%
  filter(Chromosome == 7)
view(chr7_maize_dash)
chr8_maize_dash <- maize_dash %>%
  filter(Chromosome == 8)
view(chr8_maize_dash)
chr9_maize_dash <- maize_dash %>%
  filter(Chromosome == 9)
view(chr9_maize_dash)
chr10_maize_dash <- maize_dash %>%
  filter(Chromosome == 10)
view(chr10_maize_dash)
write.table(chr1_maize_dash, file = "chr1_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE) #creates .txt file for maize sorted in descending order with dash
write.table(chr2_maize_dash, file = "chr2_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE) #Repeat for all chromosomes
write.table(chr3_maize_dash, file = "chr3_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr4_maize_dash, file = "chr4_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr5_maize_dash, file = "chr5_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr6_maize_dash, file = "chr6_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr7_maize_dash, file = "chr7_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr8_maize_dash, file = "chr8_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr9_maize_dash, file = "chr9_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr10_maize_dash, file = "chr10_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
```

## Teosinte in descending order with ? replaced with dash
```{r}
sort_desc_teosinte <- join_teosinte %>% #Create object for sorting teosinte in descending order
  arrange(as.numeric(Chromosome), desc(as.numeric(Position))) #sort teosinte in descending position
view(sort_desc_teosinte) #Check work
teosinte_dash <- data.frame(lapply(sort_desc_teosinte, gsub, pattern = "[?]", replacement = "-")) #replace ? within brackets because it is wildcard
view(teosinte_dash) #check work
chr1_teosinte_dash <- teosinte_dash %>% #creates sorted teosinte objects with dash
  filter(Chromosome == 1) #Filter by chromosome 1
view(chr1_teosinte_dash) #Check work
chr2_teosinte_dash <- teosinte_dash %>%  #Repeat for all chromosomes
  filter(Chromosome == 2)
view(chr2_teosinte_dash)
chr3_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 3)
view(chr3_teosinte_dash)
chr4_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 4)
view(chr4_teosinte_dash)
chr5_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 5)
view(chr5_teosinte_dash)
chr6_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 6)
view(chr6_teosinte_dash)
chr7_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 7)
view(chr7_teosinte_dash)
chr8_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 8)
view(chr8_teosinte_dash)
chr9_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 9)
view(chr9_teosinte_dash)
chr10_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 10)
view(chr10_teosinte_dash)
write.table(chr1_teosinte_dash, file = "chr1_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE) #creates .txt file for descending teosinte with dash
write.table(chr2_teosinte_dash, file = "chr2_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE) #Repeat for all chromosomes
write.table(chr3_teosinte_dash, file = "chr3_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr4_teosinte_dash, file = "chr4_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr5_teosinte_dash, file = "chr5_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr6_teosinte_dash, file = "chr6_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr7_teosinte_dash, file = "chr7_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr8_teosinte_dash, file = "chr8_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr9_teosinte_dash, file = "chr9_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr10_teosinte_dash, file = "chr10_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
```

## Data Visualization
## Load Packages
```{r}
library(ggplot2) #load visualization library
install.packages("ggpubr") #Package for better visualization
library(ggpubr) #load package
```
## SNP Distribution for Maize Across Chromosomes
```{r}
sort_maize_bin <- sort_maize %>% #Create binned object for sorted maize
  filter(Position != "unknown") %>% #Filter out unknown positions
  filter(Position != "multiple") %>% #Filter out multiple positions
  mutate(position_binned = cut(as.numeric(Position), 10)) %>% #Cut chromosomes into 10 bins and make own category
  filter(Chromosome != "multiple") %>% #Filter out multiple chromosomes
  filter(Chromosome != "unknown") %>% #Filter out unknown chromosomes
  mutate(Chromosome = fct_relevel(Chromosome, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))) #Use levels to put visuals in chromosome order
view(sort_maize_bin) #Check work on last column for binning
maize <- ggplot(data = sort_maize_bin) + #Name visual maize
  geom_bar(mapping = aes(x = position_binned, fill = Chromosome)) + #Create bar plot with x as binned position and fill as chromosome #
  xlab("Position (basepairs)") + #Update x label to Position (basepairs)
  ylab("SNP Distribution") + #Update y label to "SNP Distribution"
  theme(axis.text.x = element_text(angle = 90)) + #Make it pretty
  facet_grid(~as.numeric(Chromosome)) #Add facet to view chromosome as numeric
maize #Check work

```
## SNP Distribution for Teosinte Across Chromosomes
```{r}
sort_teosinte_bin <- sort_teosinte %>% #Create binned object for sorted teosinte
  filter(Position != "unknown") %>% #Filter out unknown positions
  filter(Position != "multiple") %>%  #Filter out multiple positions
  mutate(position_binned = cut(as.numeric(Position), 10)) %>% #Cut chromosomes into 10 bins and make own category
  filter(Chromosome != "multiple") %>% #Filter out multiple chromosomes
  filter(Chromosome != "unknown") %>% #Filter out unknown chromosomes
  mutate(Chromosome = fct_relevel(Chromosome, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))) #Use levels to put visuals in chromosome order
view(sort_teosinte_bin) #Check work
teosinte <- ggplot(data = sort_teosinte_bin) + #Name visual teosinte
  geom_bar(mapping = aes(x = position_binned, fill = Chromosome)) + #Create bar plot with x as binned position and fill as chromosome #
  xlab("Position (basepairs)") + #Update x label to Positions (basepairs)
  ylab("SNP Distribution") + #Update y label to SNP Distribution
  theme(axis.text.x = element_text(angle = 90)) + #Make it pretty
  facet_grid(~as.numeric(Chromosome)) #Add facet to view chromosome as numeric
teosinte #Check work
```
## Combine Maize and Teosinte SNP Distributions Across Chromosomes
```{r}
SNP_figure <- ggarrange(maize, teosinte,
                        labels = c("maize", "teosinte"),
                        ncol = 1, nrow = 2) #Create visual called snp figure that combines maize and teositne
SNP_figure
#View figure of combined maize and teosine SNP distributions per chromsome
ggsave("SNP_figure.pdf", plot = SNP_figure, width = 12, height = 15, units = "in", dpi = 300) #save figure as PDF
ggsave("SNP_figure.png", plot = SNP_figure, width = 12, height = 15, units = "in", dpi = 300) #save figure as png
#There are more SNPs in chromosome 1
```
## Add Multiple and Unknown Chromsomes with SNP Positions
```{r}
unknown_maize <- sort_maize %>% #Show unknown on each chromosome for maize
  filter(Chromosome == "unknown")
view(unknown_maize) #Check work
multiple_maize <- sort_maize %>% #Show multiple possibilities for each chromosome
  filter(Chromosome == "multiple")
view(multiple_maize) #Check work
mchr1 <- nrow(chr1_maize) #Create object for chromosome 1
mchr2 <- nrow(chr2_maize) #Repeat for all chromosomes
mchr3 <- nrow(chr3_maize)
mchr4 <- nrow(chr4_maize)
mchr5 <- nrow(chr5_maize)
mchr6 <- nrow(chr6_maize)
mchr7 <- nrow(chr7_maize)
mchr8 <- nrow(chr8_maize)
mchr9 <- nrow(chr9_maize)
mchr10 <- nrow(chr10_maize)
mchrmult <- nrow(multiple_maize)
mchrunkn <- nrow(unknown_maize)
unknown_teosinte <- sort_teosinte %>% #Show unknown on each chromosome for teosinte
  filter(Chromosome == "unknown")
view(unknown_teosinte) #Check work
multiple_teosinte <- sort_teosinte %>% #Show multiple possibiliites for each chromosome
  filter(Chromosome == "multiple")
view(multiple_teosinte) #Check work
tchr1 <- nrow(chr1_teosinte)#Create object for chromosome 1
tchr2 <- nrow(chr2_teosinte) #Repeat for all chromosomes
tchr3 <- nrow(chr3_teosinte)
tchr4 <- nrow(chr4_teosinte)
tchr5 <- nrow(chr5_teosinte)
tchr6 <- nrow(chr6_teosinte)
tchr7 <- nrow(chr7_teosinte)
tchr8 <- nrow(chr8_teosinte)
tchr9 <- nrow(chr9_teosinte)
tchr10 <- nrow(chr10_teosinte)
tchrmult <- nrow(multiple_teosinte)
tchrunkn <- nrow(unknown_teosinte)
SNP_pos_df <- data.frame(chromosome = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown"),
                         species = c("maize", "maize", "maize", "maize", "maize", "maize", "maize", "maize", "maize", "maize", "maize", "maize", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte", "teosinte"),
                         SNP_positions = c(mchr1, mchr2, mchr3, mchr4, mchr5, mchr6, mchr7, mchr8, mchr9, mchr10, mchrmult, mchrunkn, tchr1, tchr2, tchr3, tchr4, tchr5, tchr6, tchr7, tchr8, tchr9, tchr10, tchrmult, tchrunkn)) #Create a dataframe with Chromsome position, species, and SNP Position columns
view(SNP_pos_df) #Check work
SNP_pos_df <- SNP_pos_df %>%
  mutate(chromosome = fct_relevel(chromosome, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown"))) #Create levels for chromosome numbers to be in order
SNP_pos_bar <- ggplot(data = SNP_pos_df, aes(x = species, y = SNP_positions, fill = species)) + #Create bar plot with species on x, SNP position on y, and fill as species
  geom_bar(stat = "identity", position = position_dodge()) + #Make pretty
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Species") + #Relabel x axis as Species
  ylab("Number of SNP Positions") + #Relabel y axis as Number of SNP Positions
  facet_grid(~chromosome)
SNP_pos_bar #Check work
ggsave("SNP_pos_bar.pdf", plot = SNP_pos_bar, width = 10, height = 5, units = "in", dpi = 300) # print to pdf
ggsave("SNP_pos_bar.png", plot = SNP_pos_bar, width = 10, height = 5, units = "in", dpi = 300) # print to png
#SNP distrubtion is equal across maize and teosinte
```
## Zygosity Visualization
```{r}
geno2 <- read.delim("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2022/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt", header = TRUE, sep = "\t") #Pull new genotype file to mutate
geno2$homozygous <- apply(geno2, 1, function(x) length(which(x == ("A/A") | x == ("T/T") | x == ("C/C") | x == ("G/G")))) #Search for all homozygous nucleotides
view(geno2) #Check work
geno2 <- geno2 %>%
  mutate(homozygous = (homoA + homoT + homoC + homoG)) #Add columns at end for each homozygous type and total count
view(geno2) #Check work on last column
geno2$heterozygous <- apply(geno2, 1, function(x) length(which(x == ("A/C") | x == ("A/T") | x == ("A/G") | x == ("C/A") | x == ("C/T") | x == ("C/G") | x == ("G/A") | x == ("G/C") | x == ("G/T") | x == ("T/A") | x == ("T/C") | x == ("T/G")))) #Search for all heterozygous
view(geno2) #Check work on last column
geno2$missing <- apply(geno2, 1, function(x) length(which(x == ("?/?")))) #Find all missing data
view(geno2) #Check work
geno_totals <- geno2 %>% #Create new object for geno totals
  filter(Group %in% c("ZMMIL", "ZMMLR", "ZMMMR", "ZMPBA", "ZMPIL", "ZMPJA")) %>% #Filter by maize and teosinte groups
  mutate(species = fct_recode(Group,
                              "maize" = "ZMMIL",
                              "maize" = "ZMMLR",
                              "maize" = "ZMMMR",
                              "teosinte" = "ZMPBA",
                              "teosinte" = "ZMPIL",
                              "teosinte" = "ZMPJA")) #Create facet to display groups per species
view(geno_totals) #Check work
tapply(geno_totals$homozygous, geno_totals$Group, FUN=sum) #prints homozygous totals per group
tapply(geno_totals$heterozygous, geno_totals$Group, FUN=sum) #prints heterozygous total per group
tapply(geno_totals$missing, geno_totals$Group, FUN=sum) #prints missing totals per group
homozygous <- aggregate(geno_totals$homozygous, by = list(Group = geno_totals$Group), FUN = sum) #Aggregates homozygous by group
heterozygous <- aggregate(geno_totals$heterozygous, by = list(Group = geno_totals$Group), FUN = sum) #Aggregates heterozygous by group
missing <- aggregate(geno_totals$missing, by = list(Group = geno_totals$Group), FUN = sum) #Aggregates missing by group
list_zygosity_sums <- c(homozygous, heterozygous, missing) #Creates list of zygosity/missing sums
zygosity_sums <- as.data.frame(list_zygosity_sums) #Makes list a data frame
zygo_sums <- zygosity_sums %>% #Rename
  rename(homozygous = x) %>% #Orders homozygous in 2nd column
  rename(heterozygous = x.1) %>% #Orders heterozygous in 4th column
  rename(missing = x.2) #Orders missing in last column
view(zygo_sums) #Prints zygosity sums per group
zygo_sums2 <- zygo_sums %>% #Creates new object
  select(Group, homozygous, heterozygous, missing) %>% #Selects just the wanted columns
  mutate(species = fct_recode(Group,
                              "maize" = "ZMMIL",
                              "maize" = "ZMMLR",
                              "maize" = "ZMMMR",
                              "teosinte" = "ZMPBA",
                              "teosinte" = "ZMPIL",
                              "teosinte" = "ZMPJA")) #Creates facets for group
view(zygo_sums2) #Removes redundant group information
zygo_pivot <- zygo_sums2 %>% #Creates object to pivot data
  select(Group, species, homozygous, heterozygous, missing) %>% #Select wanted data
  pivot_longer(., cols = c(homozygous, heterozygous, missing), names_to = "Zygosity", values_to = "Percent") #Pivots to make columns Group, Species, Zygosity, and Percent
view(zygo_pivot) #Check work
zygo_species <- ggplot(zygo_pivot, aes(fill = Zygosity, y = Percent, x = species)) + #Make bar plot with x as species, y as percent and zyogisty for fill
  geom_bar(position="fill", stat="identity") + #Make it pretty
  xlab("Species") + #Rename x label as species
  ylab("Proportion") #Rename y label as proportion of zygosity
zygo_species #Check work of zygosity per species
ggsave("Zygo_Species.pdf", plot = zygo_species, width = 10, height = 5, units = "in", dpi = 300) # print to pdf
ggsave("Zygo_Species.png", plot = zygo_species, width = 10, height = 5, units = "in", dpi = 300) # print to png
zygo_groups <- ggplot(zygo_pivot, aes(fill = Zygosity, y = Percent, x = Group)) + #Create bar graph based on groups
  geom_bar(position="fill", stat="identity") +
  xlab("Group") +
  ylab("Proportion")
zygo_groups #Check work of zygosity per group
ggsave("Zygo_Groups.pdf", plot = zygo_groups, width = 10, height = 5, units = "in", dpi = 300) # print to pdf
ggsave("Zygo_Groups.png", plot = zygo_groups, width = 10, height = 5, units = "in", dpi = 300) # print to png
```
## Freestyle Visual, GC Content
```{r}
geno5 <- read.delim("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2022/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt", header = TRUE, sep = "\t") #Pull new genotype file to mutate
geno5$GC <- apply(geno5, 1, function(x) length(which(x == ("G/C") | x == ("C/G")| x == ("G/G") | x == ("C/C")))) #Filters anything that is certainly a G or a C
geno5$not_GC <- apply(geno5, 1, function(x) length(which(x == ("?/?") | x == ("A/A") | x == ("T/T") | x == ("A/C") | x == ("A/T") | x == ("A/G") | x == ("C/A") | x == ("C/T") | x == ("G/A") | x == ("G/T") | x == ("T/A") | x == ("T/C") | x == ("T/G")))) #Filters anything that is not for sure a G or C
view(geno5) #Check work on last column               
geno_gc <- geno5 %>%
  filter(Group %in% c("ZMMIL", "ZMMLR", "ZMMMR", "ZMPBA", "ZMPIL", "ZMPJA")) %>%
  mutate(species = fct_recode(Group,
                              "maize" = "ZMMIL",
                              "maize" = "ZMMLR",
                              "maize" = "ZMMMR",
                              "teosinte" = "ZMPBA",
                              "teosinte" = "ZMPIL",
                              "teosinte" = "ZMPJA")) #By filtering for group, we organize by species.              
view(geno_gc) #Check last column for species
gc_pivot <- geno_gc %>%
  select(Group, species, GC, not_GC) %>%
  pivot_longer(., cols = c(GC, not_GC), names_to = "GC_content", values_to = "count") #Creates column just for GC content
view(gc_pivot) #Check work
gc_group_plot <- ggplot(gc_pivot, aes(fill = GC_content, y = count, x = Group)) + #Creates bar plot where group is on the x and count is on the y
  xlab("Group") + #Relabel x axis to Group
  ylab("Proportion") + #Relabel y axis to Proportion
  geom_bar(position = "fill", stat = "identity") #Make Pretty
gc_group_plot #Check work
#There is a similar amount of GC content than non GC per group
gc_species_plot <- ggplot(gc_pivot, aes(fill = GC_content, y = count, x = species)) + #Changes x axis to species
  xlab("Species") +
  ylab("Proportion") +
  geom_bar(position = "fill", stat = "identity")
gc_species_plot #Check work
#GC content is similar among both species
ggsave("GC_Content_Groups.pdf", plot = gc_group_plot, width = 10, height = 5, units = "in", dpi = 300) # print to pdf
ggsave("GC_Content_Groups.png", plot = gc_group_plot, width = 10, height = 5, units = "in", dpi = 300) # print to png
ggsave("GC_Content_Species.pdf", plot = gc_species_plot, width = 10, height = 5, units = "in", dpi = 300) # print to pdf
ggsave("GC_Content_Species.png", plot = gc_species_plot, width = 10, height = 5, units = "in", dpi = 300) # print to png
```

## Visuals


![SNP Positions by Chromosome](\Users\chamb\Documents\R_Assignment\SNP_figure.png)


![SNP Distribution Between Species](\Users\chamb\Documents\R_Assignment\SNP_pos_bar.png)


![Zygosity by Species](\Users\chamb\Documents\R_Assignment\Zygo_Species.png)


![Zygosity by Groups](\Users\chamb\Documents\R_Assignment\Zygo_Groups.png)


![Freestyle: GC Content by Group](\Users\chamb\Documents\R_Assignment\GC_Content_Groups.png)


![Freestyle: GC Content by Species](\Users\chamb\Documents\R_Assignment\GC_Content_Species.png)


