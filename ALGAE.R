### This code produces Table 1, Table 3, Table 4, Figure 4, and Figure 6 in the manuscript: 
### Burgess SC, Turner AM, Johnston EC. Niche breadth and divergence 
# in sympatric cryptic coral species (Pocillopora spp.) 
# across habitats within reefs and among algal symbionts
# Code finalized June 2024
# Code written by Scott Burgess
# Any comments or error reporting, please contact Scott Burgess: sburgess@bio.fsu.edu

# sessionInfo()
# R version 4.4.0 (2024-04-24)
# Platform: x86_64-apple-darwin20
# Running under: macOS Sonoma 14.5

# Load libraries
library('dplyr')# v1.1.4
library('tidyr') # v1.3.1
library('glmmTMB')# v1.1.9
library('DHARMa') # v0.4.6
library('vegan') # v2.6-6.1 
library('ggplot2') # v2_3.5.1
library('forcats') # v1.0.0
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library('pairwiseAdonis') #v0.4.1



####################################################
# Import data
# ITS2 Type Profile
ITS2_Profiles <- read.csv("257_20230414T204804_DBV_20230415T005314.profiles.relative.abund_and_meta.csv") 
# Convert sample_name from dashes to underscore 
ITS2_Profiles$sample_name <- gsub("-","_",ITS2_Profiles$sample_name)

# ITS2 sequence diversity (DIV)
ITS2_post_med <- read.csv("257_20230414T204804_DBV_20230415T005314.seqs.relative.abund_and_meta.csv")
# Convert sample_name from dashes to underscore 
ITS2_post_med$sample_name <- gsub("-","_",ITS2_post_med$sample_name)

# mtORF, pbsA, and distribution meta data
d <- read.csv("Distribution data.csv")
####################################################



############ Prepare and check data ################

# Convert d$mtORF.RFLP into Pocillopora species names
unique(d$mtORF.RFLP)
Convert.names <- function(dat) {ifelse(dat$mtORF.RFLP=="Haplotype 1a_Pm" |
                                         dat$mtORF.RFLP=="Haplotype 8a",
                                       "P. meandrina",
                                       ifelse(dat$mtORF.RFLP=="Haplotype 3a" |
                                                dat$mtORF.RFLP=="Haplotype 3b" |
                                                dat$mtORF.RFLP=="Haplotype 3e" |
                                                dat$mtORF.RFLP=="Haplotype 3f" |
                                                dat$mtORF.RFLP=="Haplotype 3h" ,
                                              "P. verrucosa",
                                              ifelse(dat$mtORF.RFLP=="Haplotype 2" |
                                                       dat$mtORF.RFLP=="Haplotype 11",
                                                     "P. cf. effusa",
                                                     ifelse(dat$mtORF.RFLP=="Haplotype 5a","P. acuta",
                                                            ifelse(dat$mtORF.RFLP=="Haplotype 1a_Pgra","P. grandis",
                                                                   ifelse(dat$mtORF.RFLP=="Haplotype 10", "P. tuahiniensis","Unknown")))
                                              )
                                       ))}

d$Species <- Convert.names(d)
# Check
with(d,table(mtORF.RFLP,Species))


# Examine the _rep samples in ITS2_Profiles
ITS2_Profiles[grepl("_rep",ITS2_Profiles$sample_name),2]
ITS2_Profiles %>% filter(sample_name %in% c('21_01_01_18',
                                            '21_01_01_18_rep',
                                            '21_05_01_28',
                                            '21_05_01_28_rep',
                                            '21_05_05_01',
                                            '21_05_05_01_rep',
                                            '21_20_04_06',
                                            '21_20_04_06_rep')) %>% 
  select_if(~max(.)>0)
                         
# Examine the _rep samples in ITS2_post_med
ITS2_post_med[grepl("_rep",ITS2_post_med$sample_name),2]
ITS2_post_med %>% filter(sample_name %in% c('21_01_01_18',
                                            '21_01_01_18_rep',
                                            '21_05_01_28',
                                            '21_05_01_28_rep',
                                            '21_05_05_01',
                                            '21_05_05_01_rep',
                                            '21_20_04_06',
                                            '21_20_04_06_rep')) %>% 
  select_if(~max(.)>0)


# Add Depth.m, Site, Pocillopora species, Cladocopium species, Cladocopium clade 
# to ITS2_post_med
# Note, this will remove the "_rep" samples
# names(d)
d.cols.to.use <- c(4,5,6,12,9,10)
ITS2_post_med <- merge(ITS2_post_med, d[,d.cols.to.use], by.x='sample_name', by.y='Coral.ID')

# dim(ITS2_post_med)

# Remove samples that were not assigned a Type Profile (n=3)
tmp <- ITS2_Profiles %>% select("A1fi.A1.A1nh":"D1.D4.D6.D2.2.D1h")
ind <- which(rowSums(tmp)==0)
No_Profile <- ITS2_Profiles[ind,1]
ITS2_post_med <- (ITS2_post_med %>% filter(!(sample_name %in% No_Profile)))
# dim(ITS2_post_med)

# Add Depth.m, Site, Pocillopora species names to ITS2_Profiles
# Note, this will remove the "_rep" samples
# names(d)
d.cols.to.use <- c(4,5,6,12)
ITS2_Profiles <- merge(ITS2_Profiles, d[,d.cols.to.use], by.x='sample_name', by.y='Coral.ID')
# dim(ITS2_Profiles)

# Calculate how many Type Profiles per sample
# NOTE, this step will remove 3 samples that had no ITS2 data
ITS2_Profiles_long <- ITS2_Profiles  %>% pivot_longer("A1fi.A1.A1nh":"D1.D4.D6.D2.2.D1h",
                                                      names_to = "Type_Profile",
                                                      values_to = "Relative_abundance") %>%
  filter(Relative_abundance > 0) # Remove zero values. 
# View(ITS2_Profiles_long)


Type_Profiles_per_sample <- ITS2_Profiles_long %>%
  group_by(sample_name) %>%
  summarise(total = length(Relative_abundance))
# View(Type_Profiles_per_sample)
dim(Type_Profiles_per_sample)

# Get the dominant Type Profile for each sample
ITS2_Profiles_long_dominant <- ITS2_Profiles_long %>% 
  group_by(sample_name) %>%
  top_n(1,Relative_abundance)
dim(ITS2_Profiles_long_dominant)
# View(ITS2_Profiles_long_dominant)

# Add the dominant Type Profile to ITS2_post_med
# names(ITS2_Profiles_long_dominant)
cols.to.use <- c(1,6)
ITS2_post_med <- merge(ITS2_post_med, d[,cols.to.use], by.x='sample_name', by.y='Coral.ID')

# Add dominant Type Profiles to ITS2_post_med
# names(ITS2_Profiles_long_dominant)
ind <- c(1,6)
ITS2_post_med <- merge(ITS2_post_med,ITS2_Profiles_long_dominant[,ind], by=c('sample_name'))

# Symbiont species names for each dominant Type Profile
Profile_names <- data.frame(profiles = sort(unique(ITS2_post_med$Type_Profile)),
                            Symbiontnames = c("S. microadriaticum",
                                              "S. microadriaticum",
                                              "C. pacificum",
                                              "C. pacificum",
                                              "C. latusorum",
                                              "C. latusorum",
                                              "C. latusorum",
                                              "Cladocopium spp.",
                                              "Cladocopium spp.",
                                              "Cladocopium spp.",
                                              "C. pacificum",
                                              "C. pacificum",          
                                              "C. pacificum",
                                              "C. pacificum",                         
                                              "C. pacificum",                    
                                              "C. pacificum",        
                                              "C. pacificum",                   
                                              "C. pacificum",              
                                              "C. latusorum",   
                                              "C. latusorum",              
                                              "C. latusorum",         
                                              "C. latusorum",        
                                              "C. latusorum",          
                                              "C. latusorum",              
                                              "C. latusorum",          
                                              "C. latusorum",                
                                              "C. latusorum", 
                                              "C. latusorum",          
                                              "C. latusorum",
                                              "C. latusorum",
                                              "D. glynnii",
                                              "D. glynnii",
                                              "D. glynnii",
                                              "D. glynnii",
                                              "D. glynnii",
                                              "D. glynnii",
                                              "D. glynnii"))
# Check
# View(Profile_names)

# Add the Symbiont species name to ITS2_post_med
ITS2_post_med$Symbiontnames <- Profile_names[match(ITS2_post_med$Type_Profile,Profile_names$profiles),2]
# Check
tmp <- ITS2_post_med %>% select(Species,Type_Profile,Cladocopium_species,Symbiontnames)
# View(tmp)


# Get psbA data
psbA_data <- d %>% filter(Cladocopium_species!="")

# Add psbA_data to ITS2_Profiles_long_dominant
# names(psbA_data)
ind <- c(4,5,6,8,9,10,12)
ITS2_psbA_combined <- merge(ITS2_Profiles_long_dominant, psbA_data[,ind], by.x=c('sample_name',
                                                                                 'Depth.m',
                                                                                 'Site',
                                                                                 'Species'), by.y=c('Coral.ID',
                                                                                                    'Depth.m',
                                                                                                    'Site',
                                                                                                    'Species'),all=T)
# Add symbiont names to ITS2_psbA_combined
ITS2_psbA_combined$Symbiontnames <- Profile_names[match(ITS2_psbA_combined$Type_Profile,Profile_names$profiles),2]
ITS2_psbA_combined$Symbiontnames <- ifelse(is.na(ITS2_psbA_combined$Type_Profile), ITS2_psbA_combined$Cladocopium_species,ITS2_psbA_combined$Symbiontnames)
# View(ITS2_psbA_combined)
# Cladocopium_species is based on psbA.
# Symbiontnames is based on ITS2 profiles, for samples with ITS2. Otherwise they are Cladocopium_species (psbA)

# Checks
checks <- ITS2_psbA_combined %>% select(Species,
                                        Type_Profile,
                                        Cladocopium_species,
                                        Cladocopium_clade,
                                        Symbiontnames,
                                        Depth.m,
                                        Site)
# View(checks)

# View(ITS2_psbA_combined)
###############################################




###############################################
# Sample sizes presented in TABLE 1 and Results
# mtORF
nrow(d) # 724 colonies were identified to species
with(d,table(Depth.m,Site)) # Table 1 mtORF by Depth and Site
with(d,table(Depth.m)) # Table 1 mtORF by Depth Totals
d %>% group_by(Species) %>% summarise(n=n()) %>% mutate(freq=(n/sum(n))*100) %>% arrange(desc(n))
# P. meandrina (n = 303, 41.9%)
# P. tuahiniensis (n = 182, 25.1%)
# P. verrucosa (n = 88, 12.2%)
# P. acuta (n = 75, 10.4%)
# P. grandis (n = 56, 7.73%)
# P. cf. effusa (n = 20, 2.76%)
with(d, table(Depth.m, Species))

# ITS2 
nrow(ITS2_Profiles_long_dominant) # 377 samples with ITS2
with(ITS2_Profiles_long_dominant,table(Depth.m,Site)) # Table 1 ITS2 by Depth and Site
with(ITS2_Profiles_long_dominant,table(Depth.m)) # Table 1 ITS2 by Depth Totals


# psbA
nrow(psbA_data) # 326 samples with psbA
with(psbA_data,table(Depth.m,Site)) # Table 1 pbsA by Depth and Site
with(psbA_data,table(Depth.m)) # Table 1 pbsA by Depth Totals


# Which samples had ITS2 "OR" psbA data
nrow(ITS2_psbA_combined) # 423 had ITS2 "OR" psbA data
# 423/724
# How many samples had ITS2 "AND" psbA
tmp <- ITS2_psbA_combined %>% filter(!is.na(Type_Profile), !is.na(Cladocopium_species))
nrow(tmp) # 280 had ITS2 "AND" psbA data

ITS2_psbA_combined %>% filter(!is.na(Type_Profile)) %>% summarise(n=n()) # 377 with ITS2
ITS2_psbA_combined %>% filter(!is.na(Type_Profile),is.na(Cladocopium_species)) %>% summarise(n=n()) # 97 with ITS2 did not have psbA

ITS2_psbA_combined %>% filter(!is.na(Cladocopium_species)) %>% summarise(n=n()) # 326 with psbA
ITS2_psbA_combined %>% filter(!is.na(Cladocopium_species),is.na(Type_Profile)) %>% summarise(n=n()) # 46 with ITS2 did not have psbA

##################################################






############## Table 4 ###########################
# Remember:
# Cladocopium_species is based on psbA.
# Symbiontnames is based on ITS2 profiles, for samples with ITS2. Otherwise they are Cladocopium_species (psbA)

# Which samples hosted two symbiont species?
tmp <- checks %>% filter(!is.na(Cladocopium_species))
tmp[tmp$Cladocopium_species!=tmp$Symbiontnames,] %>% arrange(Species,Symbiontnames,Cladocopium_species)
tmp[tmp$Cladocopium_species!=tmp$Symbiontnames,] %>% arrange(Species,Symbiontnames,Cladocopium_species) %>% 
  filter(Symbiontnames %in% c('C. pacificum','C. latusorum'))


############### P. meandrina
tmp <- checks %>% filter(Species=="P. meandrina")
# Any samples with two symbiont IDs? 
tmp %>% filter(!is.na(Cladocopium_species), !(Symbiontnames %in% Cladocopium_species))
nn <- nrow(tmp); nn # this many with psbA or ITS2
# C. latusorum
tmp %>% filter(Symbiontnames=="C. latusorum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# C. pacificum
tmp %>% filter(Symbiontnames=="C. pacificum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# Cladocopium_spp.
tmp %>% filter(Symbiontnames=="Cladocopium spp.") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
tmp %>% filter(Symbiontnames=="Cladocopium_spp")
# S. microadriaticum
tmp %>% filter(Symbiontnames=="S. microadriaticum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# D. glynnii
tmp %>% filter(Symbiontnames=="D. glynnii") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4


############### P. grandis 
tmp <- checks %>% filter(Species=="P. grandis")
# Any samples with two symbiont IDs? 
tmp %>% filter(!is.na(Cladocopium_species), !(Symbiontnames %in% Cladocopium_species))
nn <- nrow(tmp); nn # this many with psbA or ITS2
# C. latusorum
tmp %>% filter(Symbiontnames=="C. latusorum" | Cladocopium_species=="C. latusorum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# C. pacificum
tmp %>% filter(Symbiontnames=="C. pacificum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# Cladocopium_spp.
tmp %>% filter(Symbiontnames=="Cladocopium spp.") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
tmp %>% filter(Symbiontnames=="Cladocopium spp.")
# S. microadriaticum
tmp %>% filter(Symbiontnames=="S. microadriaticum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# D. glynnii
tmp %>% filter(Symbiontnames=="D. glynnii") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4


############### P. cf. effusa
tmp <- checks %>% filter(Species=="P. cf. effusa")
# Any samples with two symbiont IDs? 
tmp %>% filter(!is.na(Cladocopium_species), !(Symbiontnames %in% Cladocopium_species))
nn <- nrow(tmp); nn # this many with psbA or ITS2
# C. latusorum
tmp %>% filter(Symbiontnames=="C. latusorum" | Cladocopium_species=="C. latusorum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# C. pacificum
tmp %>% filter(Symbiontnames=="C. pacificum" | Cladocopium_species=="C. pacificum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# Cladocopium_spp.
tmp %>% filter(Symbiontnames=="Cladocopium spp.") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# S. microadriaticum
tmp %>% filter(Symbiontnames=="S. microadriaticum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# D. glynnii
tmp %>% filter(Symbiontnames=="D. glynnii") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4

############### P. tuahiniensis 
tmp <- checks %>% filter(Species=="P. tuahiniensis")
# Any samples with two symbiont IDs? 
tmp %>% filter(!is.na(Cladocopium_species), !(Symbiontnames %in% Cladocopium_species))
nn <- nrow(tmp); nn # this many with psbA or ITS2
# C. latusorum
tmp %>% filter(Symbiontnames=="C. latusorum" | Cladocopium_species=="C. latusorum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
tmp %>% filter(Symbiontnames=="C. latusorum" | Cladocopium_species=="C. latusorum")
# C. pacificum
tmp %>% filter(Symbiontnames=="C. pacificum" | Cladocopium_species=="C. pacificum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# Cladocopium_spp.
tmp %>% filter(Symbiontnames=="Cladocopium spp.") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
tmp %>% filter(Symbiontnames=="Cladocopium spp.")
# S. microadriaticum
tmp %>% filter(Symbiontnames=="S. microadriaticum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# D. glynnii
tmp %>% filter(Symbiontnames=="D. glynnii") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4


############### P. verrucosa 
tmp <- checks %>% filter(Species=="P. verrucosa")
# Any samples with two symbiont IDs? 
tmp %>% filter(!is.na(Cladocopium_species), !(Symbiontnames %in% Cladocopium_species))
nn <- nrow(tmp); nn # this many with psbA or ITS2
# C. latusorum
tmp %>% filter(Symbiontnames=="C. latusorum" | Cladocopium_species=="C. latusorum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
tmp %>% filter(Symbiontnames=="C. latusorum" | Cladocopium_species=="C. latusorum")
# C. pacificum
tmp %>% filter(Symbiontnames=="C. pacificum" | Cladocopium_species=="C. pacificum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# Cladocopium_spp.
tmp %>% filter(Symbiontnames=="Cladocopium spp.") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# S. microadriaticum
tmp %>% filter(Symbiontnames=="S. microadriaticum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# D. glynnii
tmp %>% filter(Symbiontnames=="D. glynnii") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
tmp %>% filter(Symbiontnames=="D. glynnii")


############### P. acuta  
tmp <- checks %>% filter(Species=="P. acuta")
# Any samples with two symbiont IDs? 
tmp %>% filter(!is.na(Cladocopium_species), !(Symbiontnames %in% Cladocopium_species))
nn <- nrow(tmp); nn # this many with psbA or ITS2
# C. latusorum
tmp %>% filter(Symbiontnames=="C. latusorum" | Cladocopium_species=="C. latusorum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# C. pacificum
tmp %>% filter(Symbiontnames=="C. pacificum" | Cladocopium_species=="C. pacificum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# Cladocopium_spp.
tmp %>% filter(Symbiontnames=="Cladocopium spp.") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
tmp %>% filter(Symbiontnames=="Cladocopium spp.")
# S. microadriaticum
tmp %>% filter(Symbiontnames=="S. microadriaticum") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4
# D. glynnii
tmp %>% filter(Symbiontnames=="D. glynnii") %>% count(Symbiontnames) %>% summarize(freq=round((sum(n)/nn)*100,2), sum=sum(n),) # Table 4



# Check C1ag
ITS2_Profiles %>% select(Species,Site,Depth.m,C1m.C1.C1ag.C42.2.C3cg.C1b.C3cw) %>% 
  filter(C1m.C1.C1ag.C42.2.C3cg.C1b.C3cw>0)
checks %>% filter(Type_Profile == 'C1m.C1.C1ag.C42.2.C3cg.C1b.C3cw') 

ITS2_Profiles %>% select(Species,Site,Depth.m,C1.C1ah.C1ag.C42.2.C3cg.C1dj.C1b) %>% 
  filter(C1.C1ah.C1ag.C42.2.C3cg.C1dj.C1b>0)
checks %>% filter(Type_Profile == 'C1.C1ah.C1ag.C42.2.C3cg.C1dj.C1b') 

ITS2_Profiles %>% select(Species,Site,Depth.m,C1.C1ag.C1ah.C42.2.C3cg.C1b.C3cw) %>% 
  filter(C1.C1ag.C1ah.C42.2.C3cg.C1b.C3cw>0)
checks %>% filter(Type_Profile == 'C1.C1ag.C1ah.C42.2.C3cg.C1b.C3cw') 

checks %>% filter(Type_Profile == 'C1.C1ag.C1ah.C42.2.C3cg.C1b.C3cw')

###########################################################################



# ################# FIGURE 4 #####################
# set the colors for plotting
Poc_colors_to_plot <- cbind.data.frame(Species = factor(c("P. tuahiniensis",
                                            "P. meandrina",
                                            "P. verrucosa",
                                            "P. grandis",
                                            "P. cf. effusa",
                                            "P. acuta")),
                         Color = c("#D55E00",
                                   "#0072B2",
                                   "#E69F00",
                                   "#56B4E9",
                                   "#009E73",
                                   "#e63946"),
                         symbol = c(4,1,5,2,3,6))

ITS2_post_med$Species <- factor(ITS2_post_med$Species,levels=c(Poc_colors_to_plot$Species))
       
# sort(unique(ITS2_post_med$Symbiontnames))
Symbiont_colors_to_plot <- cbind.data.frame(Symbiontnames = factor(c("C. latusorum",
                                                                     "C. pacificum",
                                                                     "Cladocopium spp.",
                                                                     "D. glynnii",
                                                                     "S. microadriaticum")),
                                            Color = c("#329fba",
                                                      "#b8700d",
                                                      "#a9bf4b",
                                                      "#532c94",
                                                      "#fce80d"),
                                            symbol = c(15,16,17,9,10))

# sort(unique(ITS2_post_med$Cladocopium_clade))
# Clade_colors_to_plot <- cbind.data.frame(Cladenames = factor(c("I","II","III","IV","V","VI","VII","VIII","IX")),
#                                             Color = c("#329fba",
#                                                       "#329fba",
#                                                       "#329fba",
#                                                       "#329fba",
#                                                       "#b8700d",
#                                                       "#b8700d",
#                                                       "#b8700d",
#                                                       "#b8700d",
#                                                       "#b8700d"))
# 


# ITS2_post_med$Species <- factor(ITS2_post_med$Species,levels=Poc_colors_to_plot$Species)
# levels(ITS2_post_med$Species)

####### All species
# Create distance matrix
sym_DIV <- ITS2_post_med %>% select("A1fi":"X34947_G")
sym_matrix <- vegdist(sym_DIV, method = "bray")

# Running PCoA (using cmdscale)
mod_sym <- cmdscale(sym_matrix, eig=T)

# Extract points
data.scores = data.frame(scores(mod_sym),Species = ITS2_post_med$Species)
xy <- data.scores
# plot(xy[,1:2])
# mod_sym 
# PERMANOVA
# all_test <- adonis2(sym_DIV ~ ITS2_post_med$Species, method = "bray", permutations = 9999)
# all_test

# Pairwise PERMANOVA for all Pocillopora species/haplotypes
# pw <- pairwise.adonis2(ITS2_post_med[, -c(1:10)]~Species, data=ITS2_post_med, permutations = 9999)
# pw
# Percent of variance explained by axis
percent.var <- round(mod_sym$eig*100/sum(mod_sym$eig),1)[1:2]
############################


quartz(width=8,height=4)
# quartz(width=3.2,height=3)
par(mfrow=c(1,2),mar=c(5,4,1,1))


# Add symbols and colors to post_med_no_reps
xy$symbols <- Poc_colors_to_plot[match(ITS2_post_med$Species,Poc_colors_to_plot$Species),3]
xy$col <- Poc_colors_to_plot[match(ITS2_post_med$Species,Poc_colors_to_plot$Species),2]
xy$Symbiont_symbols <- Symbiont_colors_to_plot[match(ITS2_post_med$Symbiontnames,Symbiont_colors_to_plot$Symbiontnames),3]
xy$Symbiont_col <- Symbiont_colors_to_plot[match(ITS2_post_med$Symbiontnames,Symbiont_colors_to_plot$Symbiontnames),2]
# xy$Clade_symbols <- Clade_colors_to_plot[match(ITS2_post_med$Cladocopium_clade,Clade_colors_to_plot$Cladenames),1]
# xy$Clade_col <- Clade_colors_to_plot[match(ITS2_post_med$Cladocopium_clade,Clade_colors_to_plot$Cladenames),2]

# a)
plot(xy[,1:2],type="n",las=1, 
     xlab=eval(paste(
       "PCoA1 (",percent.var[1],"%)",sep="")),
     ylab=eval(paste(
       "PCoA2 (",percent.var[2],"%)",sep="")))
mtext("a)",side=3,adj=0)

points(x=xy[,1],y=xy[,2],
         col = xy[,which(names(xy)=="col")],
         pch=xy[,which(names(xy)=="symbols")],
         cex=1)
ordiellipse(mod_sym, ITS2_post_med$Species,
            kind = "sd",
            lty = 'dashed',
            lwd = 2,
            label=F,
            cex = 0.9,
            font=2,
            col=Poc_colors_to_plot$Color)

# Add legend
legend(-0.85,-0.15, legend=Poc_colors_to_plot$Species, pch=Poc_colors_to_plot$symbol,cex=0.9, bty="n",col=Poc_colors_to_plot$Color)
# text(c(0.08,0.12,0.18,0.15,0.29,0.57),c(-0.28,-0.17,-0.05,0.3,0.15,-0.1),labels=Poc_colors_to_plot$Species,col=Poc_colors_to_plot$Color)

# b)
plot(xy[,1:2],type="n",las=1, 
     xlab=eval(paste(
       "PCoA1 (",percent.var[1],"%)",sep="")),
     ylab=eval(paste(
       "PCoA2 (",percent.var[2],"%)",sep="")))
mtext("b)",side=3,adj=0)

points(x=xy[,1],y=xy[,2],
       col = xy[,which(names(xy)=="Symbiont_col")],
       pch=xy[,which(names(xy)=="Symbiont_symbols")],
       cex=1)
ordiellipse(mod_sym, ITS2_post_med$Symbiontnames,
            kind = "sd",
            lty = 'dashed',
            lwd = 2,
            label=F,
            cex = 0.9,
            font=2,
            col=Symbiont_colors_to_plot$Color)
# Add legend
legend(-0.85,-0.2, legend=Symbiont_colors_to_plot$Symbiontnames, pch=Symbiont_colors_to_plot$symbol,cex=0.9, bty="n",col=Symbiont_colors_to_plot$Color)
###################################################### 










################# FIGURE 6 #################

make.plot <- function(x,panel_id){
  by_depth <- x
  #Rename depths
  by_depth$Depth.m <- ifelse(by_depth$Depth.m==1,"1 m", as.character(by_depth$Depth.m))
  by_depth$Depth.m <- ifelse(by_depth$Depth.m==2,"2 m", as.character(by_depth$Depth.m))
  by_depth$Depth.m <- ifelse(by_depth$Depth.m==5,"5 m", as.character(by_depth$Depth.m))
  by_depth$Depth.m <- ifelse(by_depth$Depth.m==20,"20 m", as.character(by_depth$Depth.m))
  
  #### PCoA analysis
  sym_DIV <- by_depth %>% select("A1fi":"X34947_G")
  sym_depth <- vegdist(sym_DIV, method = "bray")
  mod_sym <- cmdscale(sym_depth,eig=T)
  
  # Extract points
  data.scores = data.frame(scores(mod_sym),Depth.m = by_depth$Depth.m, sample_name = by_depth$sample_name)
  xy <- data.scores

###########
    # Remove and re-run without the outliers
  OutLiers <- xy %>% filter(Dim2 < -0.18) %>% select(sample_name)
  by_depth <- by_depth %>% filter(!sample_name %in% OutLiers$sample_name)
  sym_DIV <- by_depth %>% select("A1fi":"X34947_G")
  sym_depth <- vegdist(sym_DIV, method = "bray")
  mod_sym <- cmdscale(sym_depth,eig=T)
  data.scores = data.frame(scores(mod_sym),Depth.m = by_depth$Depth.m, sample_name = by_depth$sample_name)
  xy <- data.scores
###########  
  # plot(xy[,1:2])
  
  # PERMANOVA
  panova <- adonis2(sym_DIV ~ by_depth$Depth.m, method = "bray", permutations = 99999)
  stats <- round(c(panova$Df[1],panova$F[1],panova$'Pr(>F)'[1]),3)
  
  # Pairwise PERMANOVA for depths
  # pw <- pairwise.adonis2(by_depth[,-c(1:10)]~Depth.m, data=by_depth, permutations = 99999)
  # pw
  
  # foo <- xy %>% filter(Dim2 > 0.2) %>% select(sample_name)
  # post_med_no_reps %>% filter(sample_name %in% foo$sample_name)
  
  ## Plot 
  # Percent of variance explained by axis
  percent.var <- round(mod_sym$eig*100/sum(mod_sym$eig),1)[1:2]
  
  plot(xy[,1:2],type="n",las=1,cex.axis=1, cex.lab=1.5,
       xlab=eval(paste(
         "PCoA1 (",percent.var[1],"%)",sep="")),
       ylab=eval(paste(
         "PCoA2 (",percent.var[2],"%)",sep="")))
  
  # Add symbols and colors
  xy$symbols <- my_depth_symbols[match(by_depth$Depth.m,my_depth_symbols$Depth),2]
  xy$col <- my_depth_symbols[match(by_depth$Depth.m,my_depth_symbols$Depth),3]

  # xy$clades <- by_depth[match(xy$sample_name,by_depth$sample_name),which(names(by_depth)=="Cladocopium_clade")]
  # text(x=xy[,1],y=xy[,2],
  #      labels=xy$clades,
  #      col=adjustcolor(xy[,which(names(xy)=="col")],alpha.f = 0.6))
  
  points(x=xy[,1],y=xy[,2],
         col = adjustcolor(xy[,which(names(xy)=="col")],alpha.f = 0.4),
         pch=xy[,which(names(xy)=="symbols")],
         cex=1.5)
  ordiellipse(mod_sym, by_depth$Depth.m, 
              kind = "sd",
              # conf = 0.95,
              lty = 'dashed', 
              lwd = 2,
              label=T, 
              cex = 0.9, 
              font=2,
              col=my_depth_symbols[match(levels(factor(by_depth$Depth.m)), my_depth_symbols$Depth),3])
  
  # Add legend
  # legend('topright', legend=paste(my_depth_symbols$Depth,c("(fringing)","(backreef)","(forereef)","(forereef)")), pch = my_depth_symbols$symbol,cex=1.1, bty="n",col=my_depth_symbols$color)
  if(stats[3]==0) mtext(side=3,adj=0,line=0.5,cex=0.8,xpd=T,
                        eval(paste("n = ",nrow(by_depth),
                                   "; F = ",stats[2],
                                   "; df = ",stats[1],
                                   "; p < 0.001",sep="")))
  
  if(stats[3]!=0) mtext(side=3,adj=0,line=0.5,cex=0.8,xpd=T,
                        eval(paste("n = ",nrow(by_depth),
                                   "; F = ",stats[2],
                                   "; df = ",stats[1],
                                   "; p = ",stats[3],sep="")))
  # mtext(side=3,adj=0,line=2.5,
  #       eval(paste(panel_id,unique(by_depth$Species))))
  # text(c(0.05,0.1,0.1,0.13,0.25,0.6),c(-0.25,-0.15,-0.05,0.3,0.15,-0.1),labels=my_coralID_symbols$Species,col=my_coralID_symbols$color)
} # end of make.plot function


# Assign symbols and colors to a each depth.
# my_depth_symbols <- cbind.data.frame(Depth=c("1m","2m","5m","20m"), 
#                                      symbol=c(15,16,17,18),
#                                      color=c("#772e25","#c44536","#197278","#283d3b"))
my_depth_symbols <- cbind.data.frame(Depth=c("2 m","5 m","20 m"), 
                                     symbol=c(16,17,18),
                                     color=c("#c44536","#197278","#283d3b"))





# Make plot
quartz(width=8,height=4)
par(mfrow=c(1,3),mar=c(5,5,5,1))

#### P. tuahiniensis by depth ####
## Prepare data
by_depth <- ITS2_post_med %>% 
  filter(Cladocopium_species=="C. pacificum" | grepl("C1d.",Type_Profile),
         !grepl("D1.",Type_Profile),
         Species=="P. tuahiniensis", Depth.m!="5")
# c('21_02_01_41','21_02_01_61') # outliers

# Checks
by_depth %>% select(Cladocopium_species,Type_Profile) %>% count(Cladocopium_species,Type_Profile)
unique(by_depth$Type_Profile)
by_depth %>% count(Cladocopium_species,Depth.m)
by_depth %>% count(Depth.m)
# Depth.m  n
# 2 21
# 20 59
by_depth %>% count(Cladocopium_clade)


make.plot(by_depth, panel_id="a) ")
mtext(side=3,adj=0,line=3,
      expression(paste("a) ", italic('Cladocopium pacificum'))))
mtext(side=3,adj=0,line=1.8,
      expression(paste(italic("- Pocillopora tuahiniensis"))))


#### P. meandrina by depth ####
## Prepare data
by_depth <- ITS2_post_med %>% 
  filter(Cladocopium_species=="C. latusorum" |
           grepl('C42a',Type_Profile)
           & !(grepl("A1",Type_Profile) & grepl("D1",Type_Profile)),
         Species=="P. meandrina",
         !Depth.m %in% c("1"))
# '21_20_02_27','21_20_05_08' outliers
# Checks
by_depth %>% select(Cladocopium_species,Type_Profile) %>% count(Cladocopium_species,Type_Profile)
unique(by_depth$Type_Profile)
by_depth %>% count(Cladocopium_species,Depth.m)
by_depth %>% count(Depth.m)
# Depth.m  n
# 2 44
# 5 85
# 20 11
by_depth %>% count(Cladocopium_clade)

make.plot(by_depth, panel_id="b)")
mtext(side=3,adj=0,line=3,
      expression(paste("b) ", italic('Cladocopium latusorum'))))
mtext(side=3,adj=0,line=1.8,
      expression(paste(italic("- Pocillopora meandrina"))))

# Add legend
par(mar=c(0,0,0,0))
plot(c(0,1),c(0,1),type="n",xlab="", ylab="",axes=F)
legend(0,0.7, 
       legend=paste(my_depth_symbols$Depth,c("back reef","fore reef","fore reef")), 
       pch = my_depth_symbols$symbol,
       cex=1.5, 
       bty="n",
       col=adjustcolor(my_depth_symbols$color,alpha.f=0.4))
####################################################################









#################### Suppl. Fig. 4. ####################
foo <- ITS2_post_med %>% 
  select(Cladocopium_species,Symbiontnames,Type_Profile,Cladocopium_clade) %>% 
  filter(Cladocopium_species %in% c("C. pacificum","C. latusorum")) 
# Use samples identified as same species by both psbA and ITS2
foo <- foo[foo$Cladocopium_species==foo$Symbiontnames,]
foo$Cladocopium_clade <- factor(foo$Cladocopium_clade,levels=c("I","II","III","IV","V","VI","VII","VIII","IX"))
# Calculate frequencies
foo <- foo %>% count(Cladocopium_clade,Type_Profile) 
# View(foo)
foo <- foo %>% group_by(Cladocopium_clade) %>% arrange(desc(n),.by_group = T)

quartz(width=6,height=4)
ggp <- foo %>% 
  mutate(Type_Profile = factor(Type_Profile,levels=Type_Profile)) %>% 
  ggplot(aes(Cladocopium_clade, Type_Profile)) + 
  geom_tile(aes(fill = n)) +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  geom_vline(xintercept = 4.5) +
  xlab("Cladocopium clade") + 
  ylab("ITS2 Type Profile") +
  theme_bw()
ggp
############################################################

#############################################################
