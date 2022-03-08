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
snp_cut <- snp_pos %>%
  select(SNP_ID, Chromosome, Position) %>%
  arrange(SNP_ID)
view(snp_cut)
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
sort_maize <- join_maize %>%
  arrange(as.numeric(Chromosome), as.numeric(Position))
# sorts by chromosome number and position, in ascending order
view(sort_maize)
chr1_maize <- sort_maize %>%
  filter(Chromosome == 1)
view(chr1_maize)
chr2_maize <- sort_maize %>%
  filter(Chromosome == 2)
chr3_maize <- sort_maize %>%
  filter(Chromosome == 3)
chr4_maize <- sort_maize %>%
  filter(Chromosome == 4)
chr5_maize <- sort_maize %>%
  filter(Chromosome == 5)
chr6_maize <- sort_maize %>%
  filter(Chromosome == 6)
chr7_maize <- sort_maize %>%
  filter(Chromosome == 7)
chr8_maize <- sort_maize %>%
  filter(Chromosome == 8)
chr9_maize <- sort_maize %>%
  filter(Chromosome == 9)
chr10_maize <- sort_maize %>%
  filter(Chromosome == 10)
view(chr10_maize)
sort_teosinte <- join_teosinte %>%
  arrange(as.numeric(Chromosome), as.numeric(Position)) #arrange teosinte object by chromosome and position
view(sort_teosinte)
chr1_teosinte <- sort_teosinte %>%
  filter(Chromosome == 1) #create sorted teosinte objects
view(chr1_teosinte)
chr2_teosinte <- sort_teosinte %>%
  filter(Chromosome == 2)
chr3_teosinte <- sort_teosinte %>%
  filter(Chromosome == 3)
chr4_teosinte <- sort_teosinte %>%
  filter(Chromosome == 4)
chr5_teosinte <- sort_teosinte %>%
  filter(Chromosome == 5)
chr6_teosinte <- sort_teosinte %>%
  filter(Chromosome == 6)
chr7_teosinte <- sort_teosinte %>%
  filter(Chromosome == 7)
chr8_teosinte <- sort_teosinte %>%
  filter(Chromosome == 8)
chr9_teosinte <- sort_teosinte %>%
  filter(Chromosome == 9)
chr10_teosinte <- sort_teosinte %>%
  filter(Chromosome == 10)
view(chr10_teosinte)
write.table(chr1_maize, file = "chr1_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE) #creates file for maize
view(chr1_maize.txt)
write.table(chr2_maize, file = "chr2_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr3_maize, file = "chr3_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr4_maize, file = "chr4_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr5_maize, file = "chr5_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr6_maize, file = "chr6_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr7_maize, file = "chr7_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr8_maize, file = "chr8_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr9_maize, file = "chr9_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr10_maize, file = "chr10_maize.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr1_teosinte, file = "chr1_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE) #creates files for teosinte
write.table(chr2_teosinte, file = "chr2_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr3_teosinte, file = "chr3_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr4_teosinte, file = "chr4_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr5_teosinte, file = "chr5_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr6_teosinte, file = "chr6_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr7_teosinte, file = "chr7_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr8_teosinte, file = "chr8_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr9_teosinte, file = "chr9_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr10_teosinte, file = "chr10_teosinte.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
sort_desc_maize <- join_maize %>%
  arrange(as.numeric(Chromosome), desc(as.numeric(Position))) #sort maize in descending position
view(sort_desc_maize)
maize_dash <- data.frame(lapply(sort_desc_maize, gsub, pattern = "[?]", replacement = "-")) #replace ? within brackets because it is wildcard
view(maize_dash)
chr1_maize_dash <- maize_dash %>% #creates sorted maize objects with dash
  filter(Chromosome == 1)
view(chr1_maize_dash)
chr2_maize_dash <- maize_dash %>%
  filter(Chromosome == 2)
chr3_maize_dash <- maize_dash %>%
  filter(Chromosome == 3)
chr4_maize_dash <- maize_dash %>%
  filter(Chromosome == 4)
chr5_maize_dash <- maize_dash %>%
  filter(Chromosome == 5)
chr6_maize_dash <- maize_dash %>%
  filter(Chromosome == 6)
chr7_maize_dash <- maize_dash %>%
  filter(Chromosome == 7)
chr8_maize_dash <- maize_dash %>%
  filter(Chromosome == 8)
chr9_maize_dash <- maize_dash %>%
  filter(Chromosome == 9)
chr10_maize_dash <- maize_dash %>%
  filter(Chromosome == 10)
view(chr10_maize_dash)
write.table(chr1_maize_dash, file = "chr1_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE) #creates file for maize dash
write.table(chr2_maize_dash, file = "chr2_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr3_maize_dash, file = "chr3_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr4_maize_dash, file = "chr4_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr5_maize_dash, file = "chr5_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr6_maize_dash, file = "chr6_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr7_maize_dash, file = "chr7_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr8_maize_dash, file = "chr8_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr9_maize_dash, file = "chr9_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr10_maize_dash, file = "chr10_maize_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
sort_desc_teosinte <- join_teosinte %>%
  arrange(as.numeric(Chromosome), desc(as.numeric(Position))) #sort teosinte in descending position
view(sort_desc_teosinte)
teosinte_dash <- data.frame(lapply(sort_desc_teosinte, gsub, pattern = "[?]", replacement = "-")) #replace ? within brackets because it is wildcard
view(teosinte_dash)
chr1_teosinte_dash <- teosinte_dash %>% #creates sorted teosinte objects with dash
  filter(Chromosome == 1)
view(chr1_teosinte_dash)
chr2_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 2)
chr3_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 3)
chr4_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 4)
chr5_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 5)
chr6_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 6)
chr7_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 7)
chr8_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 8)
chr9_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 9)
chr10_teosinte_dash <- teosinte_dash %>% 
  filter(Chromosome == 10)
view(chr10_teosinte_dash)
write.table(chr1_teosinte_dash, file = "chr1_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE) #creates file for teosinte dash
write.table(chr2_teosinte_dash, file = "chr2_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr3_teosinte_dash, file = "chr3_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr4_teosinte_dash, file = "chr4_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr5_teosinte_dash, file = "chr5_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr6_teosinte_dash, file = "chr6_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr7_teosinte_dash, file = "chr7_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr8_teosinte_dash, file = "chr8_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr9_teosinte_dash, file = "chr9_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(chr10_teosinte_dash, file = "chr10_teosinte_dash.txt", sep = "\t", row.names = FALSE, col.names = TRUE)
library(ggplot2) #load visualization library
ggplot(data = sort_maize) + geom_bar(mapping=aes(x=Chromosome))
ggplot(data = sort_teosinte) + geom_bar(mapping=aes(x=Chromosome)) #testing different variations
sort_maize_bin <- sort_maize %>%
  filter(Position != "unknown") %>%
  filter(Position != "multiple") %>%
  mutate(position_binned = cut(as.numeric(Position), 10)) %>%
  filter(Chromosome != "multiple") %>%
  filter(Chromosome != "unknown")
ggplot(data = sort_maize_bin) +
  geom_bar(mapping = aes(x = position_binned, fill = Chromosome))+
  facet_grid(~as.numeric(Chromosome)) +
  xlab("Position (basepairs)") + ylab("SNP Distribution") 

