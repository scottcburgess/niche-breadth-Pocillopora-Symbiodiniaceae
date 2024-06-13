### This code produces Figure 7, Figure 8 in the manuscript: 
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

library('dplyr')# v1.1.4
library('tidyr') # v1.3.1


# Function to calculate Levin's niche breadth for each sample
# https://www.zoology.ubc.ca/~krebs/downloads/krebs_chapter_14_2021.pdf
Levins <- function(x){
  # x is a vector of proportions across a resource type (habitat or DIV)
  # Note some formulations are (1/R) / sum(x^2), where R is the number of environments
  B <- 1 / sum(x^2) 
  (B - 1) / (length(x) - 1) # convert to scale of 0-1
  # Will give more weight to the abundant resources, compared to Shannon-Weiner
}

# Function to calculate Shannon-Weiner niche breadth for each sample
Shannon <- function(x){
  # x is a vector of proportions across a resource type (habitat or DIV)
  # won't work when any x = 0
  H <- -sum(x * log(x))
  H / log(length(x)) # convert to scale of 0-1
  # Will give more weight to the rare resources, compared to Levins
}

# Function to calculate Horn-Morisita index of niche overlap
Horn_Morisita <- function(p1,p2){
  # p1 <- c(0.5,0.5,0,0,0)
  # p2 <- c(0,0,0,0.5,0.5)
  # horn_morisita(p1,p2)
  # (2 * sum(p1*p2)) / (sum(p1^2) + sum(p2^2))
    (2 * sum(p1*p2)) / (sum(p1^2) + sum(p2^2))
}

Horn <- function(p1,p2){
  X = 1 # because we have proportion data, not abundance
  Y = X
  (sum(p1+p2) * log(p1+p2) - sum(p1*log(p1)) - sum(p2*log(p2))) / 
    (X+Y)*log(X+Y) - X*log(X) - Y*log(Y)
}   
# p1 <- c(0.5,0.5,0,0,0)
# p2 <- c(0,0,0,0.5,0.5)
# Horn(p1,p2)


# set the colors for plotting
cols <- cbind.data.frame(Species = factor(c("P. tuahiniensis",
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


######### Read in the data
# mtORF, pbsA, and distribution meta data
d <- read.csv("Distribution data.csv")

# ITS2 sequence diversity (DIV)
ITS2_post_med <- read.csv("257_20230414T204804_DBV_20230415T005314.seqs.absolute.abund_and_meta.csv")
# Convert sample_name from dashes to underscore 
ITS2_post_med$sample_name <- gsub("-","_",ITS2_post_med$sample_name)





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
d$Species <- factor(d$Species, levels=cols$Species)


# Add Depth.m, Site, Pocillopora species, Cladocopium species, 
# to ITS2_post_med
# Note, this will remove the "_rep" samples
# names(d)
d.cols.to.use <- c(4,5,6,12,9)
ITS2_post_med <- merge(ITS2_post_med, d[,d.cols.to.use], by.x='sample_name', by.y='Coral.ID')
# dim(ITS2_post_med)

# Remove samples that were not assigned a Type Profile and low depth (n=4)
tmp <- ITS2_post_med %>% select("A1fi":"X34947_G")
ITS2_post_med <- ITS2_post_med[which(rowSums(tmp)>1000),]
dim(ITS2_post_med)

# Remove sequences (columns) with all 0's
remove <- colSums(tmp)[colSums(tmp) < 1] # X535440_C
ITS2_post_med <- ITS2_post_med %>% select(!names(remove))

# Convert absolute abundance to relative abundance
tmp <- ITS2_post_med %>% select("A1fi":"X34947_G")
ITS2_post_med[,which(names(ITS2_post_med)=="A1fi"):which(names(ITS2_post_med)=="X34947_G")] <- tmp / rowSums(tmp)
# View(ITS2_post_med)

# How many ITS2 DIVs
tmp <- ITS2_post_med %>% select("A1fi":"X34947_G")
ncol(tmp) # 475 

# Calculate number of samples of each species at each site*habitat
# Prepare 'presence/absence' data for glm models
d$prob_Ptua <- ifelse(d$Species=="P. tuahiniensis",1,0)
d$prob_Pmea <- ifelse(d$Species=="P. meandrina",1,0)
d$prob_Pgra <- ifelse(d$Species=="P. grandis",1,0)
d$prob_Pver <- ifelse(d$Species=="P. verrucosa",1,0)
d$prob_Peff <- ifelse(d$Species=="P. cf. effusa",1,0)
d$prob_Pacu <- ifelse(d$Species=="P. acuta",1,0)

# # Calculate raw proportions for each species at each site*habitat
Raw_Results_Ptua <- aggregate(prob_Ptua ~ Site * Depth.m, FUN=mean, data=d)
names(Raw_Results_Ptua)[3] <- "fit"
Raw_Results_Ptua$Species <- "P. tuahiniensis"
Raw_Results_Pmea <- aggregate(prob_Pmea ~ Site * Depth.m, FUN=mean, data=d)
names(Raw_Results_Pmea)[3] <- "fit"
Raw_Results_Pmea$Species <- "P. meandrina"
Raw_Results_Pver <- aggregate(prob_Pver ~ Site * Depth.m, FUN=mean, data=d)
names(Raw_Results_Pver)[3] <- "fit"
Raw_Results_Pver$Species <- "P. verrucosa"
Raw_Results_Pgra <- aggregate(prob_Pgra ~ Site * Depth.m, FUN=mean, data=d)
names(Raw_Results_Pgra)[3] <- "fit"
Raw_Results_Pgra$Species <- "P. grandis"
Raw_Results_Peff <- aggregate(prob_Peff ~ Site * Depth.m, FUN=mean, data=d)
names(Raw_Results_Peff)[3] <- "fit"
Raw_Results_Peff$Species <- "P. cf. effusa"
Raw_Results_Pacu <- aggregate(prob_Pacu ~ Site * Depth.m, FUN=mean, data=d)
names(Raw_Results_Pacu)[3] <- "fit"
Raw_Results_Pacu$Species <- "P. acuta"
raw_props <- rbind.data.frame(Raw_Results_Ptua,
                              Raw_Results_Pmea,
                              Raw_Results_Pgra,
                              Raw_Results_Pver,
                              Raw_Results_Peff,
                              Raw_Results_Pacu)

# Calculate "relative abundance" based on the proportions of each species at each site*habitat
df_habitat <- raw_props %>% group_by(Species, Site) %>% 
  mutate(rel_abund = fit / sum(fit))

# Convert from long to wide format
df_habitat <- df_habitat %>% select(Species,Site,Depth.m,rel_abund) %>% 
  pivot_wider(names_from = Depth.m, values_from = rel_abund)
# View(df_habitat)
df_habitat <- data.frame(df_habitat)


####### Create data frame with Levin's niche breadths for habitat
Niche_habitat <- df_habitat %>% select(Species,Site)
Niche_habitat$cols <- cols[match(Niche_habitat$Species,cols$Species),2]
Niche_habitat$Levins_nb <- apply(df_habitat[,-1:-2],1,FUN=Levins)
Niche_habitat$Species <- factor(Niche_habitat$Species,levels=c("P. acuta",
                                      "P. cf. effusa",
                                      "P. tuahiniensis",
                                      "P. verrucosa",
                                      "P. grandis",
                                      "P. meandrina"))

####### Create data frame with Levin's niche breadths for symbiont species
Niche_DIV <- ITS2_post_med %>% select(Species,Depth.m,Site)
Niche_DIV$cols <- cols[match(Niche_DIV$Species,cols$Species),2]
Niche_DIV$Levins_nb <- apply(ITS2_post_med %>% select("A1fi":"X34947_G"),1,FUN=Levins)
Niche_DIV$Species <- factor(Niche_DIV$Species, levels=levels(Niche_habitat$Species))
levels(factor(Niche_DIV$Species))


# Calculate the mean niche breadth per Pocillopora species
Species_means <- Niche_habitat %>% group_by(Species) %>% 
  summarise(Levins_nb_habitat=mean(Levins_nb),
            Levins_nb_habitatSE=sd(Levins_nb)/sqrt(length(Levins_nb)))

foo <- Niche_DIV %>% group_by(Species) %>% 
  summarise(Levins_nb_DIV=mean(Levins_nb),
            Levins_nb_DIVSE=sd(Levins_nb)/sqrt(length(Levins_nb)))
Species_means <- merge(Species_means,foo,by='Species')
Species_means$cols <-  cols[match(Species_means$Species,cols$Species),2]
Species_means$pchs <-  cols[match(Species_means$Species,cols$Species),3]
# levels(factor(Species_means$Species))

# Calculate Spearman Rank Correlation 
out <- with(Species_means, cor.test(Levins_nb_habitat,Levins_nb_DIV, method='spearman',exact=T))
result <- data.frame(stat=round(out$estimate,3), p.value=round(out$p.value,3),row.names = NULL)



########## Create data frame with Horn-Morisita overlap
species.vec <- levels(factor(df_habitat$Species))
species.vec <- unique(df_habitat$Species)
site.vec <- levels(factor(df_habitat$Site))
Overlap_dat <- expand.grid(Overlap_species = species.vec, Site = site.vec)
Overlap_dat$P_tuahiniensis <- NA
Overlap_dat$P_meandrina <- NA
Overlap_dat$P_grandis <- NA
Overlap_dat$P_verrucosa <- NA
Overlap_dat$P_cf_effusa <- NA
Overlap_dat$P_acuta <- NA

# Reef habitat
for(i in 1:length(site.vec)){
  foo <- df_habitat %>% filter(Site==site.vec[i]) %>% ungroup()
  for(j in 1:length(species.vec)){
    p1 <- foo  %>% filter(Species==species.vec[j]) %>% select(X1:X20)
    p2_mat <- foo %>% select(X1:X20)
    dat <- apply(p2_mat,1, function(x) Horn_Morisita(p1=p1,p2=x))
    Overlap_dat[Overlap_dat$Site==site.vec[i],2+j] <- dat
}}





################## Make FIGURE 7 #####################  
cls <- cols[match(levels(factor(Niche_habitat$Species)),cols$Species),]
quartz(width=9,height=3)
par(mfrow=c(1,4),mar=c(6,5,2,1),oma=c(0.5,0,0,1))
with(Niche_habitat, plot(Levins_nb ~ factor(Species),
                         xlab="", ylab="", xaxt="n", yaxt="n",
                         col=cls$Color))
text(1:6,rep(-0.1,6),levels(Niche_habitat$Species),srt=45,xpd=T,cex=1.2,adj=1)
axis(side=1,at=1:6,labels=NA,las=1,cex.axis=1.2)
axis(side=2,at=seq(0,1,0.2),las=1,cex.axis=1.2)
mtext(side=2,"Niche breadth",line=3)
mtext(side=3,"a) Reef Habitat",adj=0)

# with(Species_means,points(Levins_nb_habitat ~ factor(Species),pch=19))

with(Niche_DIV, plot(Levins_nb ~ factor(Species),
                         xlab="", ylab="", xaxt="n", yaxt="n",
                         col=cls$Color))
text(1:6,rep(-0.001,6),levels(Niche_DIV$Species),srt=45,xpd=T,cex=1.2,adj=1)
axis(side=1,at=1:6,labels=NA,las=1,cex.axis=1.2)
axis(side=2,at=seq(0,1,0.005),las=1,cex.axis=1.2)
mtext(side=2,"Niche breadth",line=4)
mtext(side=3,"b) Symbionts (ITS2 DIV)",adj=0)

# with(Species_means,points(Levins_nb_DIV ~ factor(Species),pch=19))


plot(c(0,0.7),c(0.002,0.015),type="n",xlab="", ylab="", xaxt="n", yaxt="n")
with(Species_means, points(Levins_nb_habitat,Levins_nb_DIV,
                         col=cols,
                         pch=pchs,
                         lwd=2,
                         cex=2))
with(Species_means, segments(
  Levins_nb_habitat - Levins_nb_habitatSE,
  Levins_nb_DIV,
  Levins_nb_habitat + Levins_nb_habitatSE,
  Levins_nb_DIV,
  col=cols))
with(Species_means, segments(
  Levins_nb_habitat,
  Levins_nb_DIV - Levins_nb_DIVSE,
  Levins_nb_habitat,
  Levins_nb_DIV + Levins_nb_DIVSE,
  col=cols))
mtext(side=1,"Habitat Niche breadth",line=3)
axis(side=1,at=seq(0,1,0.2),las=1,cex.axis=1.2)
mtext(side=2,"Symbiont Niche breadth",line=4)
axis(side=2,at=seq(0,1,0.005),las=1,cex.axis=1.2)
mtext(side=3,"c) Niche breadth correlation",adj=0)

# Add stats
legend('topleft',
       bty="n",
       cex=1,
       legend=c(substitute(paste(rho,"=",stat),result),
                substitute(paste(p,"=",p.value),result)))


# Add legend
plot(1,1,type="n",ylab="",xlab="",axes=F)
legend(0.15,1.4, 
       legend = Species_means$Species, 
       pch = Species_means$pchs,
       cex=1.2,
       bty="n",
       col = Species_means$cols,
       xpd=T)
##################################################### 






################## Make FIGURE 8 #####################  
cls <- cols[match(levels(factor(Overlap_dat$Overlap_species)),cols$Species),]
quartz(width=3,height=7)
par(mfrow=c(6,1),mar=c(1,0,2,1),oma=c(6,5,0,1))

y <- Overlap_dat %>% select("Overlap_species","P_tuahiniensis")
y[,2] <- ifelse(y[,2]=="1",NA,y[,2])
plot(y[,1],y[,2],col=cols[cols$Species=="P. tuahiniensis",2], xlab="", ylab="", xaxt="n", yaxt="n",ylim=c(0,1))
axis(side=1,at=1:6,labels=NA,cex.axis=1.2)
axis(side=2,at=seq(0,1,0.2),las=1,cex.axis=1.2)
mtext(side=3,"a) P. tuahiniensis",adj=0,cex=0.8)

y <- Overlap_dat %>% select("Overlap_species","P_meandrina")
y[,2] <- ifelse(y[,2]=="1",NA,y[,2])
plot(y[,1],y[,2],col=cols[cols$Species=="P. meandrina",2], xlab="", ylab="", xaxt="n", yaxt="n",ylim=c(0,1))
axis(side=1,at=1:6,labels=NA,cex.axis=1.2)
axis(side=2,at=seq(0,1,0.2),las=1,cex.axis=1.2)
mtext(side=3,"b) P. meandrina",adj=0,cex=0.8)

y <- Overlap_dat %>% select("Overlap_species","P_grandis")
y[,2] <- ifelse(y[,2]=="1",NA,y[,2])
plot(y[,1],y[,2],col=cols[cols$Species=="P. grandis",2], xlab="", ylab="", xaxt="n", yaxt="n",ylim=c(0,1))
axis(side=1,at=1:6,labels=NA,cex.axis=1.2)
axis(side=2,at=seq(0,1,0.2),las=1,cex.axis=1.2)
mtext(side=3,"c) P. grandis",adj=0,cex=0.8)

y <- Overlap_dat %>% select("Overlap_species","P_verrucosa")
y[,2] <- ifelse(y[,2]=="1",NA,y[,2])
plot(y[,1],y[,2],col=cols[cols$Species=="P. verrucosa",2], xlab="", ylab="", xaxt="n", yaxt="n",ylim=c(0,1))
axis(side=1,at=1:6,labels=NA,cex.axis=1.2)
axis(side=2,at=seq(0,1,0.2),las=1,cex.axis=1.2)
mtext(side=3,"d) P. verrucosa",adj=0,cex=0.8)

y <- Overlap_dat %>% select("Overlap_species","P_cf_effusa")
y[,2] <- ifelse(y[,2]=="1",NA,y[,2])
plot(y[,1],y[,2],col=cols[cols$Species=="P. cf. effusa",2], xlab="", ylab="", xaxt="n", yaxt="n",ylim=c(0,1))
axis(side=1,at=1:6,labels=NA,cex.axis=1.2)
axis(side=2,at=seq(0,1,0.2),las=1,cex.axis=1.2)
mtext(side=3,"e) P. cf. effusa",adj=0,cex=0.8)

y <- Overlap_dat %>% select("Overlap_species","P_acuta")
y[,2] <- ifelse(y[,2]=="1",NA,y[,2])
plot(y[,1],y[,2],col=cols[cols$Species=="P. acuta",2], xlab="", ylab="", xaxt="n", yaxt="n",ylim=c(0,1))
axis(side=1,at=1:6,labels=NA,cex.axis=1.2)
axis(side=2,at=seq(0,1,0.2),las=1,cex.axis=1.2)
mtext(side=3,"f) P. acuta",adj=0,cex=0.8)

text(1:6,rep(-0.2,6),levels(Overlap_dat$Overlap_species),srt=45,xpd=NA,cex=1.2,adj=1)

mtext(side=2,"Reef habitat niche overlap (Horn-Morisita index)",outer=T,line=3)
