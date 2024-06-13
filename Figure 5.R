### This code produces Figure 5 in the manuscript: 
### Burgess SC, Turner AM, Johnston EC. Niche breadth and divergence 
# in sympatric cryptic coral species (Pocillopora spp.) 
# across habitats within reefs and among algal symbionts
# Code finalized June 2024
# Code adapted by Erika Johnston from code originally written by: 
# M.Marzonie and M. Nitschke and found here: 
# https://github.com/magenamarzonie/CoralSeaSymbiont/blob/main/1_Analysis.Rmd
# Any comments or error reporting, please contact Scott Burgess: sburgess@bio.fsu.edu

library(ggsci)
library(tidyverse)
library(dplyr)
library(forcats)
library(reshape2)
library(stringr)
library(tidyr)
library(tibble)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biostrings")
library(Biostrings)
library(phangorn)
library(ape)
library(ggplot2)
library(ggtree)
library(patchwork)
library(bioseq)
library(kmer)
library(GUniFrac)
library(seqinr)
library(vegan)
library(corrplot)
library(ggrepel)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ggmsa")
library(ggmsa)
library(dendextend)
library(usedist)

# Load symbiont data
seqs <- read.csv("257_20230414T204804_DBV_20230415T005314.seqs.absolute.abund_and_meta.csv")

seqs$sample_name <- gsub("-","_",seqs$sample_name)

# Load metadata
meta <- read.csv("Distribution data.csv") %>%
  filter(Cladocopium_clade != "") %>%
  mutate(Habitat=case_when(Depth.m=="1"~"Fringing (1m)",
                           Depth.m=="2"~"Back reef (2m)",
                           Depth.m=="5"~"Fore reef (5m)",
                           Depth.m=="20"~"Fore reef (20m)")) %>%
  select(-Depth.m) %>%
  mutate(Poc_species=case_when(mtORF.RFLP=="Haplotype 10"~"P. tuahiniensis",
                               mtORF.RFLP=="Haplotype 11"~"P. cf. effusa",
                               mtORF.RFLP=="Haplotype 2"~"P. cf. effusa",
                               mtORF.RFLP=="Haplotype 1a_Pm"~"P. meandrina",
                               mtORF.RFLP=="Haplotype 8a"~"P. meandrina",
                               mtORF.RFLP=="Haplotype 1a_Pgra"~"P. grandis",
                               mtORF.RFLP=="Haplotype 3a"~"P. verrucosa",
                               mtORF.RFLP=="Haplotype 3b"~"P. verrucosa",
                               mtORF.RFLP=="Haplotype 3e"~"P. verrucosa",
                               mtORF.RFLP=="Haplotype 3f"~"P. verrucosa",
                               mtORF.RFLP=="Haplotype 3h"~"P. verrucosa",
                               mtORF.RFLP=="Haplotype 5a"~"P. acuta")) %>%
  select(-mtORF.RFLP,-Trip,-Date.collected.YYYY.MM.DD,-Preservative,-Site,-Box.ID,-Plate.Number) %>%
  dplyr::rename(sample_name=Coral.ID)

# Filter seq data for specific samples
samples <- meta %>%
  select(sample_name)

seqs <- filter(seqs, sample_name %in% samples$sample_name)


##### Load Custom Functions
read_fasta_df <- function (file = "") {
  fasta <- readLines(file)
  ind <- grep(">", fasta)
  s <- data.frame(ind = ind, from = ind + 1, to = c((ind - 
                                                       1)[-1], length(fasta)))
  seqs <- rep(NA, length(ind))
  for (i in 1:length(ind)) {
    seqs[i] <- paste(fasta[s$from[i]:s$to[i]], collapse = "")
  }
  tib <- tibble(label = gsub(">", "", fasta[ind]), sequence = seqs)
  return(tib)
}

write_fasta_df <- function (data, filename) 
{
  fastaLines = c()
  for (rowNum in 1:nrow(data)) {
    fastaLines = c(fastaLines, as.character(paste(">", 
                                                  data[rowNum, "label"], sep = "")))
    fastaLines = c(fastaLines, as.character(data[rowNum, 
                                                 "sequence"]))
  }
  fileConn <- file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

dna_to_DNAbin <- function (dna){
  DNAbin <- as_DNAbin(dna)
  names(DNAbin) <- names(dna)
  return(DNAbin)
}
dna_to_DNAStringset <- function(x) 
{
  bioseq:::check_dna(x)
  DNAstr <- DNAStringSet(paste(x))
  names(DNAstr) <- names(x)
  return(DNAstr)
}

DNAStringSet_to_dna <- function(x){
  x_dna <- as_dna(paste(x))
  names(x_dna) <- names(x)
  res <- tibble(label = names(x), sequence = x_dna)
  return(res)
}

# Convert DNAstringset to DNAbin
DNAStringSet_to_DNAbin <- function(DNAStringSet){
  DNAbin <- as.DNAbin(DNAStringSet)
  return(DNAbin)
}

# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2315-y
palette <- c("A" = "#46ff2d", 
             "G" = "#ffae01", 
             "C" = "#f24641", 
             "T" = "#4294fa", 
             "K" = "#8b4816",
             "M" = "#83831f",
             "R" = "#ffff81",
             "S" = "#ff9d80",
             "Y" = "#e381f2",
             "W" = "#80fff2",
             "V" = "#fde4b8",
             "B" = "#f9c1bf",
             "H" = "#c0d9f9",
             "D" = "#c7ffba",
             "U" = "#8989fb",
             "N" = "black", 
             "-" = "white",
             "+" = "White")


pal_df <- data.frame(names = names(palette), col = palette)


# Convert to long format
seqs_long <- seqs %>%
  filter(!is.na(sample_name)) %>%
  select(sample_name, "A1fi":"X34947_G") %>%
  pivot_longer("A1fi":"X34947_G") %>%
  filter(value > 0)# Remove zero values


seqs_long <- left_join(seqs_long, meta, by="sample_name") 

# Q. Are we working with the post-med seqs according to the metadata in seqs?
san_check <- seqs_long %>%
  group_by(sample_name) %>%
  summarise(total = sum(value)) #A. yes

# Create a list of samples to keep that didn't fail to sequence
keepers_ss <- san_check %>%
  filter(total > 1000)

non_keep <- san_check %>% 
  filter(total < 1000)

#we filter out 4 samples
# Filter out the failed samples
seqs_long <- seqs_long %>%
  filter(sample_name %in% keepers_ss$sample_name) %>%
  group_by(sample_name) %>%
  mutate(value_rel = value/sum(value)) %>% # Convert to relative abundance
  ungroup() %>%
  mutate(name = as.factor(name)) # Make sample names a factor


# Create a random palette for each sequence
n <- length(levels(seqs_long$name))
seqs_pal = rainbow(n, s=.6, v=.9)[sample(1:n,n, replace = FALSE)]
names(seqs_pal) <- levels(seqs_long$name)

# Read in the profile data
profiles_raw <- read.csv("257_20230414T204804_DBV_20230415T005314.profiles.absolute.abund_and_meta.csv") %>%
  select(sample_name, `A1fi.A1.A1nh`:`D1.D4.D6.D2.2.D1h`)

# Convert sample_name from dashes to underscore 
profiles_raw$sample_name <- gsub("-","_",profiles_raw$sample_name)

#Convert to long format 
profiles_long <- profiles_raw %>%
  pivot_longer(`A1fi.A1.A1nh`:`D1.D4.D6.D2.2.D1h`) %>% # Convert it to long format
  mutate(name = paste0("p_", name)) %>% # Add a p_ to the beginning of each profile (Some profiles are single sequence profiles and clash with the Sequence names)
  filter(sample_name %in% seqs_long$sample_name) %>% # Remove samples that dont appear in the Sequence dataframe
  group_by(sample_name) %>%
  mutate(value = as.numeric(value)) %>%
  filter(value > 0) %>% # Remove 0 abundance profiles
  mutate(sample_name = as.factor(sample_name),
         name = as.factor(name)) %>% 
  ungroup() %>%
  left_join(., meta, by="sample_name") # Add in metadata

# What is the total number of profile-related sequences in each sample?
profiles_sum <- profiles_long %>%
  group_by(sample_name) %>%
  summarise(total = sum(value))

# How many sequences in each sample are not part of a profile?
residual <- left_join(profiles_sum, san_check, by = "sample_name") %>%
  mutate(residual = total.y - total.x) %>%
  select(sample_name, value = residual) %>%
  mutate(name = "non-profile sequences") %>%
  left_join(., meta)

# Combine the profiles and non-profile sequences
profile_data <- rbind(profiles_long, residual) %>%
  group_by(sample_name) %>%
  mutate(value_rel = value/sum(value)) # convert to relative abundance - in that sample 

# Create palette for profiles
#Order ITS2 type profiles
profile_data$name = factor(profile_data$name,levels=c("p_C1d",
                                  "p_C1d.C42.2.C1.C1k.C1b",
                                  "p_C1d.C42.2.C1.C1k.C1b.C3cg",
                                  "p_C1d.C42g.C42.2.C1.C1b.C42a",
                                  "p_C1d.C1.C42.2.C3.C1b.C3cg.C45c.C115k",
                                  "p_C1d.C1.C42.2.C1b.C3cg.C45c.C115k.C1au",
                                  "p_C1d.C42.2.C1.C3cg.C1b.C3cw.C45c.C115k",
                                  "p_C1d.C1.C42.2.C3.C1b.C3cg.C45c.C115k.C1au.C41p",
                                  "p_C3.C115d",
                                  "p_C1m.C1.C1ag.C42.2.C3cg.C1b.C3cw",
                                  "p_C1.C1ah.C1ag.C42.2.C3cg.C1dj.C1b",
                                  "p_C1.C1ag.C1ah.C42.2.C3cg.C1b.C3cw",
                                  "p_C1.C42.2.C1b.C3.C42g.C42a.C1az",
                                  "p_C42a.C42.2.C1b.C1.C42b.C1au",
                                  "p_C42a.C42.2.C1.C1b.C1j.C1au.C3",
                                  "p_C42a.C1.C42.2.C1b.C1j.C1au.C3.C115l",
                                  "p_C42a.C42.2.C1.C42b.C1b.C1az.C1au.C3",
                                  "p_C42a.C42.2.C1.C1b.C42ad.C42ac.C1az.C1au.C42g",
                                  "p_C42.2.C1.C1b.C3.C1au.C41p.C115k",
                                  "p_C42g.C1.C42.2.C42a.C1b.C1az.C1au.C3",
                                  "p_C42.2.C1.C42g.C42a.C1b.C42b.C3.C1au",
                                  "p_C1.C42.2.C1b.C1az.C42g.C1au.C3.C42az",
                                  "p_C42.2.C1.C42a.C1b.C1au.C1az.C3.C115d",
                                  "p_C42.2.C1.C42a.C42g.C1b.C42b.C1az.C1au",
                                  "p_C1.C42.2.C42g.C42a.C1b.C1au.C1az.C3.C42h",
                                  "p_C42.2.C1.C1b.C1au.C3.C1az.C115k.C41p.C115d",
                                  "p_C116.C116i.C116q",
                                  "p_D1.D2d.D1ky.D1kx",
                                  "p_D1.D2d.D1aa.D1z.D1y",
                                  "p_D1.D2d.D1aa.D1z.D1hx",
                                  "p_D1.D2d.D1aa.D1z.D1y.D1x",
                                  "p_D1.D6.D4.D2.2.D1h.D2f.D2c",
                                  "p_D1.D6.D4.D2d.D2.2.D1h.D2f",
                                  "p_A1.A1fi",
                                  "p_A1fi.A1.A1nh"))

#Set ITS2 profile colors
profile_pal <- c("p_C1d"="#fad8b9",
                       "p_C1d.C42.2.C1.C1k.C1b"="#edb179",
                       "p_C1d.C42.2.C1.C1k.C1b.C3cg"="#e69f5e",
                       "p_C1d.C42g.C42.2.C1.C1b.C42a"="#e6ad5e",
                       "p_C1d.C1.C42.2.C3.C1b.C3cg.C45c.C115k"="#de9f47",
                       "p_C1d.C1.C42.2.C1b.C3cg.C45c.C115k.C1au"="#d18f32",
                       "p_C1d.C42.2.C1.C3cg.C1b.C3cw.C45c.C115k"="#c27d1d",
                       "p_C1d.C1.C42.2.C3.C1b.C3cg.C45c.C115k.C1au.C41p"="#b8700d",
                       "p_C3.C115d"="#f5abab",
                       "p_C1m.C1.C1ag.C42.2.C3cg.C1b.C3cw"="#f77c7c",
                       "p_C1.C1ah.C1ag.C42.2.C3cg.C1dj.C1b"="#db5a5a",
                       "p_C1.C1ag.C1ah.C42.2.C3cg.C1b.C3cw"="#a62314",
                       "p_C1.C42.2.C1b.C3.C42g.C42a.C1az"="#cafaef",
                       "p_C42a.C42.2.C1b.C1.C42b.C1au"="#9feddb",
                       "p_C42a.C42.2.C1.C1b.C1j.C1au.C3"="#7cebd1",
                       "p_C42a.C1.C42.2.C1b.C1j.C1au.C3.C115l"="#5ee6c6",
                       "p_C42a.C42.2.C1.C42b.C1b.C1az.C1au.C3"="#3dd9b4",
                       "p_C42a.C42.2.C1.C1b.C42ad.C42ac.C1az.C1au.C42g"="#1fbf99",
                       "p_C42.2.C1.C1b.C3.C1au.C41p.C115k"="#5ad1d1",
                       "p_C42g.C1.C42.2.C42a.C1b.C1az.C1au.C3"="#5eb2b8",
                       "p_C42.2.C1.C42g.C42a.C1b.C42b.C3.C1au"="#3aa8b0",
                       "p_C1.C42.2.C1b.C1az.C42g.C1au.C3.C42az"="#329fba",
                       "p_C42.2.C1.C42a.C1b.C1au.C1az.C3.C115d"="#288da6",
                       "p_C42.2.C1.C42a.C42g.C1b.C42b.C1az.C1au"="#11728a",
                       "p_C1.C42.2.C42g.C42a.C1b.C1au.C1az.C3.C42h"="#086c85",
                       "p_C42.2.C1.C1b.C1au.C3.C1az.C115k.C41p.C115d"="#02586e",
                       "p_C116.C116i.C116q"="#a9bf4b",
                       "p_D1.D2d.D1ky.D1kx"="#754cba",
                       "p_D1.D2d.D1aa.D1z.D1y"="#532c94",
                       "p_D1.D2d.D1aa.D1z.D1hx"="#371b63",
                       "p_D1.D2d.D1aa.D1z.D1y.D1x"="#1f1038",
                       "p_D1.D6.D4.D2.2.D1h.D2f.D2c"="#fc9ddc",
                       "p_D1.D6.D4.D2d.D2.2.D1h.D2f"="#fc7ed2",
                       "p_A1.A1fi"="#f5e753",
                       "p_A1fi.A1.A1nh"="#fce80d")


# Merge the palettes and replace the non-profile sequences with grey
all_pal <- c(seqs_pal, profile_pal)
all_pal['non-profile sequences'] <- "#808080" 

# Join profiles and sequence data together into single dataframe and add more metadata
all_data <- rbind(seqs_long, profile_data)


##1.1 Basic Stats
# How many samples per species?
all_data %>%
  distinct(sample_name, Poc_species) %>%
  group_by(Poc_species) %>% 
  summarise(total_samples = n())

# P. acuta                    6
# P. cf. effusa              10
# P. grandis                 21
# P. meandrina              119
# P. tuahiniensis            85
# P. verrucosa               38

# Total number of sequences (for whole library)? 
study_total <- all_data %>% 
  filter(!(str_detect(name, "p_")),
         name != "non-profile sequences") %>% 
  summarise(total_seqs = sum(value)) %>%
  pull(total_seqs)
#8,005,025 total 	

#Total number of sequences (per species)? 
all_data %>%
  filter(!(str_detect(name, "p_")),
         name != "non-profile sequences") %>% 
  group_by(Poc_species) %>%
  summarise(total_seqs = sum(value))

# P. acuta            204111
# P. cf. effusa       444480
# P. grandis          484921
# P. meandrina       3277082
# P. tuahiniensis    2723260
# P. verrucosa        871171

##1.2 Sequencing Depth
# Average per sample sequencing depth? (per coral species)

# Average per sample across all spp
all_data %>%
  group_by(sample_name) %>% 
  dplyr::summarise(per_sample = sum(value)) %>%
  dplyr::summarise(mean_all = mean(per_sample))
# Across all:  57,384	

#average per sample per coral host species 
all_data %>% 
  group_by(Poc_species, sample_name) %>% 
  summarise(per_sample = sum(value)) %>%
  summarise(mean_all = mean(per_sample), n = n())

# P. acuta          68037      6
# P. cf. effusa     88896     10
# P. grandis        46183.    21
# P. meandrina      55077.   119
# P. tuahiniensis   64077.    85
# P. verrucosa      45851.    38	


##1.3 Number of sequences per symbiont genus 
# Total number of Cladocopium 
all_data %>%
  filter(!(str_detect(name, "p_")), name != "non-profile sequences") %>%
  filter(str_sub(name, 1, 1) == "C" | str_detect(name, "_C")) %>% 
  summarise(sum = sum(value))

(7692283 / study_total) * 100 # 96.09318 Clado

#Total number of Symbiodinium
all_data %>%
  filter(!(str_detect(name, "p_")), name != "non-profile sequences") %>%
  filter(str_sub(name, 1, 1) == "A" | str_detect(name, "_A")) %>% 
  summarise(sum = sum(value)) 

(44150 / study_total) * 100 # 0.5515286 Symbio


# Total number of Durusdinium sequences
all_data %>%
  filter(!(str_detect(name, "p_")), name != "non-profile sequences") %>%
  filter(str_sub(name, 1, 1) == "D" | str_detect(name, "_D")) %>% 
  summarise(sum = sum(value)) 

(268394 / study_total) * 100 # 3.352819 Durusdinium


##1.4 Total number of Type Profiles 
#total type profiles across Pocilloporidae 
all_data %>%
  filter(str_detect(name, "p_")) %>%    #profiles start with p_
  group_by(name) %>%
  dplyr:: count() %>%
  dplyr:: arrange(desc(n)) %>%
  ungroup() %>%
  mutate(prop = n/sum(n)) %>%
  mutate(cumulative_sum = cumsum(prop)) %>%
  print(n=34)

# profile                                                n    prop cumulative_sum
# p_C42g.C1.C42.2.C42a.C1b.C1az.C1au.C3              76 0.258            0.258
# p_C1d.C1.C42.2.C3.C1b.C3cg.C45c.C115k.C1au.C41p    62 0.210            0.468
# p_C1.C42.2.C42g.C42a.C1b.C1au.C1az.C3.C42h         43 0.146            0.614
# p_C1d.C42.2.C1.C1k.C1b.C3cg                        19 0.0644           0.678
# p_C1.C1ag.C1ah.C42.2.C3cg.C1b.C3cw                 13 0.0441           0.722
# p_C1d.C42.2.C1.C3cg.C1b.C3cw.C45c.C115k             9 0.0305           0.753
# p_C42.2.C1.C1b.C3.C1au.C41p.C115k                   8 0.0271           0.780
# p_C1m.C1.C1ag.C42.2.C3cg.C1b.C3cw                   6 0.0203           0.8  
# p_C42a.C1.C42.2.C1b.C1j.C1au.C3.C115l               5 0.0169           0.817
# p_C42.2.C1.C42a.C42g.C1b.C42b.C1az.C1au             5 0.0169           0.834
# p_D1.D2d.D1aa.D1z.D1hx                              5 0.0169           0.851
# p_C42.2.C1.C42g.C42a.C1b.C42b.C3.C1au               4 0.0136           0.864
# p_D1.D2d.D1aa.D1z.D1y.D1x                           4 0.0136           0.878
# p_C1d.C1.C42.2.C1b.C3cg.C45c.C115k.C1au             3 0.0102           0.888
# p_C42a.C42.2.C1.C1b.C1j.C1au.C3                     3 0.0102           0.898
# p_C1.C42.2.C1b.C1az.C42g.C1au.C3.C42az              3 0.0102           0.908
# p_D1.D6.D4.D2d.D2.2.D1h.D2f                         3 0.0102           0.919
# p_C1d.C42.2.C1.C1k.C1b                              2 0.00678          0.925
# p_C1d.C42g.C42.2.C1.C1b.C42a                        2 0.00678          0.932
# p_C1.C1ah.C1ag.C42.2.C3cg.C1dj.C1b                  2 0.00678          0.939
# p_C42a.C42.2.C1.C42b.C1b.C1az.C1au.C3               2 0.00678          0.946
# p_C42.2.C1.C1b.C1au.C3.C1az.C115k.C41p.C115d        2 0.00678          0.953
# p_D1.D2d.D1ky.D1kx                                  2 0.00678          0.959
# p_A1fi.A1.A1nh                                      2 0.00678          0.966
# p_C1d                                               1 0.00339          0.969
# p_C1d.C1.C42.2.C3.C1b.C3cg.C45c.C115k               1 0.00339          0.973
# p_C3.C115d                                          1 0.00339          0.976
# p_C42a.C42.2.C1b.C1.C42b.C1au                       1 0.00339          0.980
# p_C42a.C42.2.C1.C1b.C42ad.C42ac.C1az.C1au.C42g      1 0.00339          0.983
# p_C42.2.C1.C42a.C1b.C1au.C1az.C3.C115d              1 0.00339          0.986
# p_C116.C116i.C116q                                  1 0.00339          0.990
# p_D1.D2d.D1aa.D1z.D1y                               1 0.00339          0.993
# p_D1.D6.D4.D2.2.D1h.D2f.D2c                         1 0.00339          0.997
# p_A1.A1fi                                           1 0.00339          1.00    

#total number of profiles, in order of highest abundance across P. acuta samples 
all_data %>%
  filter(Poc_species == "P. acuta") %>% 
  filter(str_detect(name, "p_")) %>%    #profiles start with p_
  group_by(name) %>%
  dplyr:: count() %>%
  dplyr:: arrange(desc(n)) %>%
  ungroup() %>%
  mutate(prop = n/sum(n)) %>%
  mutate(cumulative_sum = cumsum(prop))

# p_D1.D2d.D1aa.D1z.D1y.D1x     3 0.5            0.5  
# p_D1.D2d.D1ky.D1kx            1 0.167          0.667
# p_D1.D2d.D1aa.D1z.D1y         1 0.167          0.833
# p_D1.D2d.D1aa.D1z.D1hx        1 0.167          1

#total number of profiles, in order of highest abundance across P. cf. effusa samples 
all_data %>%
  filter(Poc_species == "P. cf. effusa") %>% 
  filter(str_detect(name, "p_")) %>%
  group_by(name) %>%
  dplyr:: count() %>%
  dplyr:: arrange(desc(n)) %>%
  ungroup() %>%
  mutate(prop = n/sum(n)) %>%
  mutate(cumulative_sum = cumsum(prop))

# p_C42.2.C1.C1b.C3.C1au.C41p.C115k                   8 0.727           0.727
# p_C1d.C1.C42.2.C3.C1b.C3cg.C45c.C115k.C1au.C41p     2 0.182           0.909
# p_D1.D6.D4.D2.2.D1h.D2f.D2c                         1 0.0909          1 

#total number of profiles, in order of highest abundance across P. grandis samples 
all_data %>%
  filter(Poc_species == "P. grandis") %>% 
  filter(str_detect(name, "p_")) %>%    #profiles start with p_
  group_by(name) %>%
  dplyr:: count() %>%
  dplyr:: arrange(desc(n)) %>%
  ungroup() %>%
  mutate(prop = n/sum(n)) %>%
  mutate(cumulative_sum = cumsum(prop))

# p_C42g.C1.C42.2.C42a.C1b.C1az.C1au.C3              9  0.36           0.36
# p_C42a.C1.C42.2.C1b.C1j.C1au.C3.C115l              4  0.16           0.52
# p_C1.C42.2.C42g.C42a.C1b.C1au.C1az.C3.C42h         4  0.16           0.68
# p_C42a.C42.2.C1.C1b.C1j.C1au.C3                    2  0.08           0.76
# p_A1fi.A1.A1nh                                     2  0.08           0.84
# p_C42a.C42.2.C1b.C1.C42b.C1au                      1  0.04           0.88
# p_C42a.C42.2.C1.C1b.C42ad.C42ac.C1az.C1au.C42g     1  0.04           0.92
# p_D1.D2d.D1aa.D1z.D1y.D1x                          1  0.04           0.96
# p_A1.A1fi                                          1  0.04           1    

#total number of profiles, in order of highest abundance across P. meandrina samples 
all_data %>%
  filter(Poc_species == "P. meandrina") %>% 
  filter(str_detect(name, "p_")) %>%    #profiles start with p_
  group_by(name) %>%
  dplyr:: count() %>%
  dplyr:: arrange(desc(n)) %>%
  ungroup() %>%
  mutate(prop = n/sum(n)) %>%
  mutate(cumulative_sum = cumsum(prop))

# p_C42g.C1.C42.2.C42a.C1b.C1az.C1au.C3           67 0.540            0.540
# p_C1.C42.2.C42g.C42a.C1b.C1au.C1az.C3.C42h      37 0.298            0.839
# p_C42.2.C1.C42a.C42g.C1b.C42b.C1az.C1au          5 0.0403           0.879
# p_C1.C42.2.C1b.C1az.C42g.C1au.C3.C42az           3 0.0242           0.903
# p_D1.D2d.D1aa.D1z.D1hx                           3 0.0242           0.927
# p_C42.2.C1.C1b.C1au.C3.C1az.C115k.C41p.C115d     2 0.0161           0.944
# p_C1d                                            1 0.00806          0.952
# p_C42a.C42.2.C1.C1b.C1j.C1au.C3                  1 0.00806          0.960
# p_C42a.C1.C42.2.C1b.C1j.C1au.C3.C115l            1 0.00806          0.968
# p_C42a.C42.2.C1.C42b.C1b.C1az.C1au.C3            1 0.00806          0.976
# p_C42.2.C1.C42g.C42a.C1b.C42b.C3.C1au            1 0.00806          0.984
# p_C42.2.C1.C42a.C1b.C1au.C1az.C3.C115d           1 0.00806          0.992
# p_D1.D2d.D1ky.D1kx                               1 0.00806          1    

#total number of profiles, in order of highest abundance across P. tuahiniensis samples 
all_data %>%
  filter(Poc_species == "P. tuahiniensis") %>% 
  filter(str_detect(name, "p_")) %>%    #profiles start with p_
  group_by(name) %>%
  dplyr:: count() %>%
  dplyr:: arrange(desc(n)) %>%
  ungroup() %>%
  mutate(prop = n/sum(n)) %>%
  mutate(cumulative_sum = cumsum(prop))

# p_C1d.C1.C42.2.C3.C1b.C3cg.C45c.C115k.C1au.C41p    54 0.607           0.607
# p_C1d.C42.2.C1.C1k.C1b.C3cg                        15 0.169           0.775
# p_C1d.C42.2.C1.C3cg.C1b.C3cw.C45c.C115k             9 0.101           0.876
# p_C42.2.C1.C42g.C42a.C1b.C42b.C3.C1au               3 0.0337          0.910
# p_D1.D6.D4.D2d.D2.2.D1h.D2f                         3 0.0337          0.944
# p_C1d.C42g.C42.2.C1.C1b.C42a                        1 0.0112          0.955
# p_C1d.C1.C42.2.C3.C1b.C3cg.C45c.C115k               1 0.0112          0.966
# p_C1d.C1.C42.2.C1b.C3cg.C45c.C115k.C1au             1 0.0112          0.978
# p_C1.C42.2.C42g.C42a.C1b.C1au.C1az.C3.C42h          1 0.0112          0.989
# p_C116.C116i.C116q                                  1 0.0112          1  

#total number of profiles, in order of highest abundance across P. verrucosa samples 
all_data %>%
  filter(Poc_species == "P. verrucosa") %>% 
  filter(str_detect(name, "p_")) %>%    #profiles start with p_
  group_by(name) %>%
  dplyr:: count() %>%
  dplyr:: arrange(desc(n)) %>%
  ungroup() %>%
  mutate(prop = n/sum(n)) %>%
  mutate(cumulative_sum = cumsum(prop))

# p_C1.C1ag.C1ah.C42.2.C3cg.C1b.C3cw                 13 0.325          0.325
# p_C1d.C1.C42.2.C3.C1b.C3cg.C45c.C115k.C1au.C41p     6 0.15           0.475
# p_C1m.C1.C1ag.C42.2.C3cg.C1b.C3cw                   6 0.15           0.625
# p_C1d.C42.2.C1.C1k.C1b.C3cg                         4 0.1            0.725
# p_C1d.C42.2.C1.C1k.C1b                              2 0.05           0.775
# p_C1d.C1.C42.2.C1b.C3cg.C45c.C115k.C1au             2 0.05           0.825
# p_C1.C1ah.C1ag.C42.2.C3cg.C1dj.C1b                  2 0.05           0.875
# p_C1d.C42g.C42.2.C1.C1b.C42a                        1 0.025          0.9  
# p_C3.C115d                                          1 0.025          0.925
# p_C42a.C42.2.C1.C42b.C1b.C1az.C1au.C3               1 0.025          0.95 
# p_C1.C42.2.C42g.C42a.C1b.C1au.C1az.C3.C42h          1 0.025          0.975
# p_D1.D2d.D1aa.D1z.D1hx                              1 0.025          1   

##1.5 Type profile proportions
#What is the proportion of samples that had 1 and 2 type profiles?

all_data %>% 
  filter(str_detect(name, "p_")) %>%
  group_by(sample_name) %>% 
  summarise(n = n()) %>% 
  filter(n == 1)     #263  samples have 1 type profile
#filter(n == 2)     #16 samples have 2 type profiles


# 2 UPGMA Stats
#filter out non-profile sequences
seq_data <- all_data %>% 
  filter(!str_detect(name, "non")) %>% 
  filter(!str_detect(name, "p_"))

#creating an object for each species 
p_acuta_seqs <- seq_data %>% filter(Poc_species == "P. acuta")
p_effusa_seqs <- seq_data %>% filter(Poc_species == "P. cf. effusa")
p_grandis_seqs <- seq_data %>% filter(Poc_species == "P. grandis")
p_meandrina_seqs <- seq_data %>% filter(Poc_species == "P. meandrina")
p_tuahiniensis_seqs <- seq_data %>% filter(Poc_species == "P. tuahiniensis")
p_verrucosa_seqs <- seq_data %>% filter(Poc_species == "P. verrucosa")

#colour palette
# Merge the palettes and replace the non-profile sequences with grey
all_pal <- c(seqs_pal, profile_pal)
all_pal['non-profile sequences'] <- "#808080" 



## 2.1 Poci UPGMA Tree
library(bioseq)
library(kmer)
library(GUniFrac)
library(poppr)
library(ggtree)

####### Tree with Cladocopium, Durusdinium, and Symbiodinium #####

#read in file 
fasta_poci <- read_fasta_df("257_20230414T204804_DBV_20230415T005314.seqs.fasta") %>%
  filter(label %in% seq_data$name) %>%   #only keeping DNA seqs that appear in seqs_long subset 
  #filter(str_detect(label, "D")) %>%
  deframe() %>%
  as_dna();fasta_poci

#creating the tree
kdist_poci <- fasta_poci %>%
  dna_to_DNAbin() %>%
  kdistance(k = 7, residues = "DNA", method = "edgar") %>%
  as.matrix()

k_tree_poci <- kdist_poci %>% phangorn::upgma()

seqs_wide_poci <- seq_data %>%
  dplyr::select(sample_name, name, value) %>%
  #filter(str_detect(name, "D")) %>%
  filter(!str_detect(name, "X")) %>%
  pivot_wider(names_from = name, values_from = value, values_fill = 0) %>%
  tibble::column_to_rownames(var = "sample_name")

k_unidist_poci <- GUniFrac(seqs_wide_poci, k_tree_poci)   #GUniFrac calculates all the distances 
k_unidist_poci <- k_unidist_poci$unifracs

du_poci <- k_unidist_poci[, , "d_0.5"]    # GUniFrac with alpha 0.5 
dist_poci <- as.dist(du_poci, diag = FALSE)

# Cluster the samples
hclust_samps_poci <- upgma(du_poci)

# Make the sample tree with tip labels
tree_poci <- ggtree(hclust_samps_poci, size = 0.2) +
  theme(aspect.ratio = 1) +
  geom_tiplab(size=0, offset=0.1); tree_poci # Set size to 1 to see sample name

# Get a sample order from ggtree
poci_sample_order <- tree_poci$data %>% filter(isTip == "TRUE") %>%
  arrange(y) %>%
  pull(label)

pocimat <- meta %>% filter(sample_name %in% tree_poci$data$label) %>%
  select(sample_name, Poc_species, Habitat, Cladocopium_clade) %>%
  tibble::column_to_rownames(var = "sample_name")


poci_tree_mat <- gheatmap(tree_poci, pocimat, colnames = F, width=0.5, offset = -0.475)


poci_tree_mat <- poci_tree_mat + scale_fill_manual(values = c("lightgrey", "#424242", "darkgrey", "white",
                                                              "#94f7f4","#0d85bd","#12704f","#0d077a","#ab0308",
                                                              "#e63946","#009E73","#56B4E9","#0072B2","#D55E00","#E69F00",
                                                              "#f5cc9a","#e0b310","#f2555d","#e03f16",
                                                              "grey"))+
  layout_dendrogram()+ 
  theme(aspect.ratio = 0.3, legend.position = "right")


# Start plotting the composition data
plot_df_poci <- all_data %>%
  mutate(Vial = fct_relevel(sample_name, poci_sample_order))

theme_set(theme_bw())

## Fig 5A
bar_uni_poci <- 
  ggplot(plot_df_poci, aes(Vial, value_rel)) +
  geom_bar(stat = "identity", aes(fill = name, colour = name)) +
  theme(aspect.ratio = 0.25, legend.position = "bottom", axis.text.y=element_blank(), # change to legend.position = "bottom" to plot legend
        legend.key.spacing.y = unit(0.01, "cm"), # unmark to plot legend
        axis.ticks.y = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.ticks = element_blank()) +
  scale_fill_manual(values = all_pal, breaks = levels(profile_data$name)) +
  scale_colour_manual(values = all_pal, breaks = levels(profile_data$name)) +
  geom_hline(yintercept = 1, linewidth = 1) +
  guides(fill=guide_legend(ncol=2)); bar_uni_poci 

## Plot together
poci_tree_mat / bar_uni_poci
