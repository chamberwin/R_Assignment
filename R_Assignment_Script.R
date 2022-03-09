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
install.packages("ggpubr") #Package for better visualization
library(ggpubr)
sort_maize_bin <- sort_maize %>%
  filter(Position != "unknown") %>%
  filter(Position != "multiple") %>%
  mutate(position_binned = cut(as.numeric(Position), 10)) %>%
  filter(Chromosome != "multiple") %>%
  filter(Chromosome != "unknown") %>%
  mutate(Chromosome = fct_relevel(Chromosome, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")))
view(sort_maize_bin)
maize <- ggplot(data = sort_maize_bin) +
  geom_bar(mapping = aes(x = position_binned, fill = Chromosome)) +
  xlab("Position(basepairs)") +
  ylab("SNP Distribution") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(~as.numeric(Chromosome))
sort_teosinte_bin <- sort_teosinte %>%
  filter(Position != "unknown") %>%
  filter(Position != "multiple") %>%
  mutate(position_binned = cut(as.numeric(Position), 10)) %>%
  filter(Chromosome != "multiple") %>%
  filter(Chromosome != "unknown") %>%
  mutate(Chromosome = fct_relevel(Chromosome, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")))
teosinte <- ggplot(data = sort_teosinte_bin) +
  geom_bar(mapping = aes(x = position_binned, fill = Chromosome)) +
  xlab("Position(basepairs)") +
  ylab("SNP Distribution") +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(~as.numeric(Chromosome))
SNP_figure <- ggarrange(maize, teosinte,
                        labels = c("maize", "teosinte"),
                        ncol = 1, nrow = 2)
SNP_figure #View figure
ggsave("SNP_figure.pdf", plot = SNP_figure, width = 12, height = 15, units = "in", dpi = 300) #save figure as PDF
ggsave("SNP_figure.png", plot = SNP_figure, width = 12, height = 15, units = "in", dpi = 300)
unknown_maize <- sort_maize %>% #Show unknown on each chromosome
  filter(Chromosome == "unknown")
view(unknown_maize)
multiple_maize <- sort_maize %>% #Show multiple possibilities for each chromosome
  filter(Chromosome == "multiple")
view(multiple_maize)
mchr1 <- nrow(chr1_maize)
mchr2 <- nrow(chr2_maize)
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
unknown_teosinte <- sort_teosinte %>%
  filter(Chromosome == "unknown")
view(unknown_teosinte)
multiple_teosinte <- sort_teosinte %>%
  filter(Chromosome == "multiple")
view(multiple_teosinte)
tchr1 <- nrow(chr1_teosinte)
tchr2 <- nrow(chr2_teosinte)
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
                         SNP_positions = c(mchr1, mchr2, mchr3, mchr4, mchr5, mchr6, mchr7, mchr8, mchr9, mchr10, mchrmult, mchrunkn, tchr1, tchr2, tchr3, tchr4, tchr5, tchr6, tchr7, tchr8, tchr9, tchr10, tchrmult, tchrunkn))
view(SNP_pos_df)
SNP_pos_df <- SNP_pos_df %>%
  mutate(chromosome = fct_relevel(chromosome, c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "multiple", "unknown")))
SNP_pos_bar <- ggplot(data = SNP_pos_df, aes(x = species, y = SNP_positions, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Species") +
  ylab("Number of SNP positions") +
  facet_grid(~chromosome)
SNP_pos_bar
ggsave("SNP_pos_bar.pdf", plot = SNP_pos_bar, width = 10, height = 5, units = "in", dpi = 300) # print to pdf
ggsave("SNP_pos_bar.png", plot = SNP_pos_bar, width = 10, height = 5, units = "in", dpi = 300) # print to png
geno2 <- read.delim("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2022/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt", header = TRUE, sep = "\t") #Pull new genotype file to mutate
geno2$homozygous <- apply(geno2, 1, function(x) length(which(x == ("A/A") | x == ("T/T") | x == ("C/C") | x == ("G/G")))) #Search for all homozygous nucleotides
view(geno2)
geno2 <- geno2 %>%
  mutate(homozygous = (homoA + homoT + homoC + homoG)) #Add columns at end for each homozygous type and total count
view(geno2)
geno2$heterozygous <- apply(geno2, 1, function(x) length(which(x == ("A/C") | x == ("A/T") | x == ("A/G") | x == ("C/A") | x == ("C/T") | x == ("C/G") | x == ("G/A") | x == ("G/C") | x == ("G/T") | x == ("T/A") | x == ("T/C") | x == ("T/G"))))
view(geno2)
geno2$missing <- apply(geno2, 1, function(x) length(which(x == ("?/?")))) #Find all missing data
view(geno2)
geno_totals <- geno2 %>%
  filter(Group %in% c("ZMMIL", "ZMMLR", "ZMMMR", "ZMPBA", "ZMPIL", "ZMPJA")) %>%
  mutate(species = fct_recode(Group,
                              "maize" = "ZMMIL",
                              "maize" = "ZMMLR",
                              "maize" = "ZMMMR",
                              "teosinte" = "ZMPBA",
                              "teosinte" = "ZMPIL",
                              "teosinte" = "ZMPJA"))
view(geno_totals)
tapply(geno_totals$homozygous, geno_totals$Group, FUN=sum) #prints homozygous totals per group
tapply(geno_totals$hetezygous, geno_totals$Group, FUN=sum) #prints heterozygous total per group
tapply(geno_totals$missing, geno_totals$Group, FUN=sum) #prints missing totals per group
homozygous <- aggregate(geno_totals$homozygous, by = list(Group = geno_totals$Group), FUN = sum)
heterozygous <- aggregate(geno_totals$hetezygous, by = list(Group = geno_totals$Group), FUN = sum)
missing <- aggregate(geno_totals$missing, by = list(Group = geno_totals$Group), FUN = sum)
list_zygosity_sums <- c(homozygous, heterozygous, missing)
zygosity_sums <- as.data.frame(list_zygosity_sums)
zygo_sums <- zygosity_sums %>%
  rename(homozygous = x) %>%
  rename(heterozygous = x.1) %>%
  rename(missing = x.2)
view(zygo_sums)
zygo_sums2 <- zygo_sums %>%
  select(Group, homozygous, heterozygous, missing) %>%
  mutate(species = fct_recode(Group,
                              "maize" = "ZMMIL",
                              "maize" = "ZMMLR",
                              "maize" = "ZMMMR",
                              "teosinte" = "ZMPBA",
                              "teosinte" = "ZMPIL",
                              "teosinte" = "ZMPJA"))
view(zygo_sums2)
zygo_pivot <- zygo_sums2 %>%
  select(Group, species, homozygous, heterozygous, missing) %>%
  pivot_longer(., cols = c(homozygous, heterozygous, missing), names_to = "Zygosity", values_to = "Percent")
view(zygo_pivot)
zygo_species <- ggplot(zygo_pivot, aes(fill = Zygosity, y = Percent, x = species)) +
  geom_bar(position="fill", stat="identity") +
  xlab("Species") +
  ylab("Proportion")
ggsave("Zygo_Species.pdf", plot = zygo_species, width = 10, height = 5, units = "in", dpi = 300) # print to pdf
ggsave("Zygo_Species.png", plot = zygo_species, width = 10, height = 5, units = "in", dpi = 300) # print to png
zygo_groups <- ggplot(zygo_pivot, aes(fill = Zygosity, y = Percent, x = Group)) +
  geom_bar(position="fill", stat="identity") +
  xlab("Group") +
  ylab("Proportion")
zygo_groups
ggsave("Zygo_Groups.pdf", plot = zygo_groups, width = 10, height = 5, units = "in", dpi = 300) # print to pdf
ggsave("Zygo_Groups.png", plot = zygo_groups, width = 10, height = 5, units = "in", dpi = 300) # print to png
geno5 <- read.delim("https://raw.githubusercontent.com/EEOB-BioData/BCB546-Spring2022/main/assignments/UNIX_Assignment/fang_et_al_genotypes.txt", header = TRUE, sep = "\t") #Pull new genotype file to mutate
geno5$GC <- apply(geno5, 1, function(x) length(which(x == ("G/C") | x == ("C/G"))))
geno5$not_GC <- apply(geno5, 1, function(x) length(which(x == ("?/?") | x == ("A/A") | x == ("T/T") | x == ("C/C") | x == ("G/G") | x == ("A/C") | x == ("A/T") | x == ("A/G") | x == ("C/A") | x == ("C/T") | x == ("G/A") | x == ("G/T") | x == ("T/A") | x == ("T/C") | x == ("T/G"))))
view(geno5)                
geno_gc <- geno5 %>%
  filter(Group %in% c("ZMMIL", "ZMMLR", "ZMMMR", "ZMPBA", "ZMPIL", "ZMPJA")) %>%
  mutate(species = fct_recode(Group,
                              "maize" = "ZMMIL",
                              "maize" = "ZMMLR",
                              "maize" = "ZMMMR",
                              "teosinte" = "ZMPBA",
                              "teosinte" = "ZMPIL",
                              "teosinte" = "ZMPJA"))                
gc_pivot <- geno_gc %>%
  select(Group, species, GC, not_GC) %>%
  pivot_longer(., cols = c(GC, not_GC), names_to = "GC_content", values_to = "count")
view(gc_pivot)
gc_group_plot <- ggplot(gc_pivot, aes(fill = GC_content, y = count, x = Group)) +
  xlab("Group") +
  ylab("Proportion") +
  geom_bar(position = "fill", stat = "identity")
gc_group_plot
gc_species_plot <- ggplot(gc_pivot, aes(fill = GC_content, y = count, x = species)) +
  xlab("Species") +
  ylab("Proportion") +
  geom_bar(position = "fill", stat = "identity")
gc_species_plot
ggsave("GC_Content_Groups.pdf", plot = gc_group_plot, width = 10, height = 5, units = "in", dpi = 300) # print to pdf
ggsave("GC_Content_Groups.png", plot = gc_group_plot, width = 10, height = 5, units = "in", dpi = 300) # print to png
ggsave("GC_Content_Species.pdf", plot = gc_species_plot, width = 10, height = 5, units = "in", dpi = 300) # print to pdf
ggsave("GC_Content_Species.png", plot = gc_species_plot, width = 10, height = 5, units = "in", dpi = 300) # print to png
