library(tidyverse)
library(dplyr)
install.packages("janitor")
install.packages("data.table")
library(janitor)
library(data.table)
geno <- read.delim("https://raw.githubusercontent.com/chamberwin/BCB546-Spring2022/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt", header = TRUE, sep = "\t")
snp_pos <- read.delim("https://raw.githubusercontent.com/chamberwin/BCB546-Spring2022/main/assignments/UNIX_Assignment/snp_position.txt", header = TRUE, sep = "\t") 
ncol(geno)
dim(geno)
length(geno)
object.size(geno)
summary(geno)
nrow(snp_pos)
ncol(snp_pos)
dim(snp_pos)
length(snp_pos)
object.size(snp_pos)
summary(snp_pos)
## Maize Join
geno_maize <- geno %>%
  filter(Group %in% c("ZMMIL", "ZMMLR", "ZMMMR", "Sample_ID")) %>%
  # filter by maize groups only, including sample ID
  t() %>%
  # transpose - note it's added automatic column names (V1, V2, etc. in new header)
  as.data.frame() %>%
  # make sure it's a dataframe
  setDT(keep.rownames = TRUE) %>%
  #gives column 1 a placeholder name as well
  row_to_names(row_number = 1) %>%
  # sets first row as column names
  rename(SNP_ID = Sample_ID) # rename first column as SNP_ID so we can join on like column names
view(geno_maize) #double check joined by SNP_ID
join_maize <- snp_cut %>%
  inner_join(geno_maize, by = "SNP_ID") # joins both files - inner_join makes sure we take only SNP_IDs in both files
view(join_maize) #double check it joined correctly
##Teosinte Join
geno_teosinte <- geno %>%
  filter(Group %in% c("ZMPBA", "ZMPIL", "ZMPJA")) %>%
  # filter by maize groups only, including sample ID
  t() %>%
  # transpose - note it's added automatic column names (V1, V2, etc. in new header)
  as.data.frame() %>%
  # make sure it's a dataframe
  setDT(keep.rownames = TRUE) %>%
  #gives column 1 a placeholder name as well
  row_to_names(row_number = 1) %>%
  # sets first row as column names
  rename(SNP_ID = Sample_ID)
# rename first column as SNP_ID so we can join on like column names
view(geno_teosinte) #double check it is orded by SNP_ID
join_teosinte <- snp_cut %>%
  inner_join(geno_teosinte, by = "SNP_ID")
# joins both files - inner_join makes sure we take only SNP_IDs in both files
view(join_teosinte) #double check it joined correctly
