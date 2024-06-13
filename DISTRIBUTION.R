### This code produces Table 3, Figure 2, Figure 3, Figure S1 in the manuscript: 
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


# Function to fit a 2-way (categorical) glmmTMB model and extract results from LR ratio test
# and fitted values and 95% CI for plotting. 
# On the link scale
# For Site * Depth.m
fit.models <- function(data,y.name,family){
  # For testing:
  # data=d
  # y.name="prob_Pmea"
  # family="binomial"
  foo <- data[,which(names(data) %in% c(y.name,"Site","Depth.m"))]
  names(foo)[which(names(foo)==y.name)] <- "y"
  foo$Depth.m <- factor(foo$Depth.m)
  foo$Site <- factor(foo$Site)
  
  m1 <- glmmTMB(y ~ Site * Depth.m, family=family, data=foo)
  m2 <- glmmTMB(y ~ Site + Depth.m, family=family, data=foo)
  m3 <- glmmTMB(y ~ Site, family=family, data=foo)
  m4 <- glmmTMB(y ~ Depth.m, family=family, data=foo)
  
  predictions <- expand.grid(Site=unique(foo$Site), Depth.m=unique(foo$Depth.m)) 
  p <- predict(m1,predictions,re.form=NA,se.fit=T,type="link")
  predictions <- cbind.data.frame(predictions,p)
  predictions$lwr <- predictions$fit - 2*predictions$se.fit
  predictions$upr <- predictions$fit + 2*predictions$se.fit
  predictions$fit <- plogis(predictions$fit) 
  predictions$lwr <- plogis(predictions$lwr) 
  predictions$upr <- plogis(predictions$upr) 
  
  LRtable <- data.frame(model=c("Interaction","Additive effect of Depth.m","Additive effect of Site"),
                        ChiSq = c(round(na.omit(anova(m1,m2)$Chisq),3),
                                  round(na.omit(anova(m2,m3)$Chisq),3),
                                  round(na.omit(anova(m2,m4)$Chisq),3)),
                        df = c(na.omit(anova(m1,m2)$'Chi Df'),
                               na.omit(anova(m2,m3)$'Chi Df'),
                               na.omit(anova(m2,m4)$'Chi Df')),
                        p.val = c(round(na.omit(anova(m1,m2)$'Pr(>Chisq)'),3),
                                  round(na.omit(anova(m2,m3)$'Pr(>Chisq)'),3),
                                  round(na.omit(anova(m2,m4)$'Pr(>Chisq)'),3)))
  
  return(list(LRtable=LRtable,predictions=predictions))
}



# Function to fit a 1-way (categorical) glmmTMB model and extract results from LR ratio test
# and fitted values and 95% CI for plotting. 
# On the link scale
# For Depth.m only
fit.models.pooled <- function(data,y.name,family){
  foo <- data[,which(names(data) %in% c(y.name,"Site","Depth.m"))]
  names(foo)[which(names(foo)==y.name)] <- "y"
  foo$Depth.m <- factor(foo$Depth.m)
  foo$Site <- factor(foo$Site)
  
  m1 <- glmmTMB(y ~ Depth.m, family=family, data=foo)
  m2 <- glmmTMB(y ~ 1, family=family, data=foo)

  predictions <- expand.grid(Depth.m=unique(foo$Depth.m)) 
  p <- predict(m1,predictions,re.form=NA,se.fit=T,type="link")
  predictions <- cbind.data.frame(predictions,p)
  predictions$lwr <- predictions$fit - 2*predictions$se.fit
  predictions$upr <- predictions$fit + 2*predictions$se.fit
  predictions$fit <- plogis(predictions$fit) 
  predictions$lwr <- plogis(predictions$lwr) 
  predictions$upr <- plogis(predictions$upr) 
  
  LRtable <- data.frame(model=c("Effect of Depth.m"),
                        ChiSq = c(round(na.omit(anova(m1,m2)$Chisq),3)),
                        df = c(na.omit(anova(m1,m2)$'Chi Df')),
                        p.val = c(round(na.omit(anova(m1,m2)$'Pr(>Chisq)'),3)))
  
  return(list(LRtable=LRtable,predictions=predictions))
}





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
                                   "#e63946"))





# Read in the data
d <- read.csv("Distribution data.csv")



############ Prepare and check data ################
# Convert mtORF.RFLP into Pocillopora species names
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

# Prepare 'presence/absence' data for glm models
d$prob_Ptua <- ifelse(d$Species=="P. tuahiniensis",1,0)
d$prob_Pmea <- ifelse(d$Species=="P. meandrina",1,0)
d$prob_Pgra <- ifelse(d$Species=="P. grandis",1,0)
d$prob_Pver <- ifelse(d$Species=="P. verrucosa",1,0)
d$prob_Peff <- ifelse(d$Species=="P. cf. effusa",1,0)
d$prob_Pacu <- ifelse(d$Species=="P. acuta",1,0)

# Check structure of data
with(d, table(Depth.m,Site,Species)) # Can't estimate Depth.m * Site * Species

# How many samples from each species?
d %>% group_by(Species) %>% summarise(n=n()) %>% mutate(freq=(n/sum(n))*100) %>% arrange(desc(n))
nrow(d)
# Species             n  freq
# P. meandrina      303 41.9 
# P. tuahiniensis   182 25.1 
# P. verrucosa       88 12.2 
# P. acuta           75 10.4 
# P. grandis         56  7.73
# P. cf. effusa      20  2.76
#################################################




#### Fit models, get summaries, convert to data scale for plotting #######
# P.tuahiniensis
# No colonies at 1m, so remove that depth, add NA to $predictions
Results_Ptua <- fit.models.pooled(data=d[d$Depth.m!="1",],y.name="prob_Ptua",family="binomial")
Results_Ptua$predictions <- rbind.data.frame(Results_Ptua$predictions, 
                                             data.frame(Depth.m=rep("1",1),
                                          fit=rep(NA,1),
                                          se.fit=rep(NA,1),
                                          lwr=rep(NA,1),
                                          upr=rep(NA,1)))

# P.meandrina
Results_Pmea <- fit.models.pooled(data=d,y.name="prob_Pmea",family="binomial")

# P.grandis
Results_Pgra <- fit.models.pooled(data=d,y.name="prob_Pgra",family="binomial")

# P.verrucosa
Results_Pver <- fit.models.pooled(data=d,y.name="prob_Pver",family="binomial") 

# P. effusa
# No colonies at 1m, so remove that depth, add NA to $predictions
Results_Peff <- fit.models.pooled(data=d[d$Depth.m!="1",],y.name="prob_Peff",family="binomial")
Results_Peff$predictions <- rbind.data.frame(Results_Peff$predictions, 
                                             data.frame(Depth.m=rep("1",1),
                                                        fit=rep(NA,1),
                                                        se.fit=rep(NA,1),
                                                        lwr=rep(NA,1),
                                                        upr=rep(NA,1)))

# P.acuta
# No colonies at 5m, 10, or 20m so remove that depth, add 0's to $predictions
# d %>% count(Species,Depth.m)
Results_Pacu <- fit.models.pooled(data=d[d$Depth.m %in% c("1","2"),],y.name="prob_Pacu",family="binomial")
Results_Pacu$predictions <- rbind.data.frame(Results_Pacu$predictions, 
                                             data.frame(Depth.m=c("5","10","20"),
                                                        fit=rep(NA,3),
                                                        se.fit=rep(NA,3),
                                                        lwr=rep(NA,3),
                                                        upr=rep(NA,3)))

# Add a column for Species ID
Results_Ptua$predictions$Species <- cols$Species[1]
Results_Pmea$predictions$Species <- cols$Species[2]
Results_Pver$predictions$Species <- cols$Species[3]
Results_Pgra$predictions$Species <- cols$Species[4]
Results_Peff$predictions$Species <- cols$Species[5]
Results_Pacu$predictions$Species <- cols$Species[6]

# Combine predictions for each species for plotting
pred <- rbind.data.frame(Results_Ptua$predictions,
                         Results_Pmea$predictions,
                         Results_Pver$predictions,
                         Results_Pgra$predictions,
                         Results_Peff$predictions,
                         Results_Pacu$predictions)

# Make sure Depth.m is in correct order
pred$Depth.m <- factor(pred$Depth.m,levels=c("1","2","5","10","20"))



# # Calculate raw proportions for each site, for plotting
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
raw_props$fit <- ifelse(raw_props$fit==0,NA,raw_props$fit)
################################################################


########### Table 3 #############################
# (remove P.acuta because only one colony in back reef)
Table3 <- rbind.data.frame(Results_Ptua$LRtable,
                         Results_Pmea$LRtable,
                         Results_Pver$LRtable,
                         Results_Pgra$LRtable,
                         Results_Peff$LRtable)
Table3
####################################################





##################### FIGURE 2 #####################
xlabels=c("1m fringing","2m back reef","5m fore reef","10m fore reef","20m fore reef")

quartz(width=2.5,height=5)
par(mfrow=c(7,1),mar=c(1,0.5,0,0),oma=c(5,5,1,1))

ylims <- c(0.9,1,0.3,0.2,0.2,0.6)
bys <- c(0.3,0.3,0.1,0.05,0.05,0.2)

for(i in 1:7){
    if(i!=7){
      y <- pred %>% filter(Species==cols$Species[i])
      y <- y[order(y$Depth.m),]
      x <- barplot(y$fit, yaxt="n", col=cols$Color[i], ylim=c(0,ylims[i]),border=NA)
      axis(side=2,at=seq(0,ylims[i],bys[i]),las=1,cex.axis=1.2)
      axis(side=2,at=seq(0,ylims[i],bys[i]),labels=NA)
      segments(x,y$lwr,x,y$upr,col="black",lwd=1.5)
      axis(side=1,at=x,pos=-0.05*ylims[i],labels=NA)
    }
    if(i==7){
      y_mat <- pred %>% select(Depth.m,Species,fit) %>% spread(Depth.m,fit)
      y_mat <- y_mat[match(cols$Species,y_mat$Species),] 
      y_mat[is.na(y_mat)] <- 0
      x <- barplot(as.matrix(y_mat[,-1]),col=cols$Color,yaxt="n",xaxt="n",ylim=c(0,1.05),border=NA)
      axis(side=2,at=seq(0,1,0.2),las=1,cex.axis=1.2)
      axis(side=2,at=seq(0,1,0.2),labels=NA)
      axis(side=1,at=x,pos=-0.05,labels=NA)
      text(x,rep(-0.3,length(xlabels)),xlabels,srt=45,xpd=NA,adj=1,cex=1.2)
    }
}
mtext("Relative abundance (proportion of Pocillopora colonies)",side=2,outer=T,line=3,cex=1)
############################################################# 







######################### FIGURE 3 #########################
### Distance-based Redundancy Analysis (db-RDA).
## See, for example:
## https://stackoverflow.com/questions/51715281/dbrda-in-r-how-to-with-abundance-data-and-missing-values-for-environmental-data#51767138
## https://archetypalecology.wordpress.com/2018/02/21/distance-based-redundancy-analysis-db-rda-in-r/

y_mat <- raw_props %>% select(Depth.m,Species,Site,fit) %>% spread(Species,fit)
y_mat$Depth.m <- paste0(y_mat$Depth.m,"m")
# y_mat <- y_mat[match(cols$Species,y_mat$Species),] 
y_mat[is.na(y_mat)] <- 0
prop_mat <- y_mat[,c(-1,-2)]

####### Model 1 - Constrained by Depth 
m1 <- dbrda(prop_mat ~ Depth.m, distance="bray",sqrt.dist=F,data=y_mat)
summary(m1) # Proportion Explained: dbRDA1 = 37.38%, dbRDA2 = 32.00%

# Perform test of each factor using marginal effect
# which is the  unique effect of each variable conditional on the presence of the other variable in the model
# anova.result <- anova.cca(m1,by="margin",permutations=999999) # Result in paper. Note that the p-values will change slightly each time this line is run, more so when permutations are reduced (quicker, but less accurate).
# anova.result

sppscores(m1) <- prop_mat # add species scores

# Extract Site scores (weighted sums of species scores)
xy <- summary(m1)$sites[,1:2]

# Extract species scores
haps <- summary(m1)$species[,1:2]

# Add symbols and colors to post_med_no_reps
my_depth_symbols <- data.frame(Depth=c("1m","2m","5m","10m","20m"),
                               pch=c(0,1,2,5,6),
                               color=c("#DBC539","#49BA67","#379FAA","#455EBA","#333872"))

col.vec <- my_depth_symbols[match(y_mat$Depth.m,my_depth_symbols$Depth),3]


## Plot 
# quartz(width=4,height=4)
quartz(width=3.2,height=3)
par(mfrow=c(1,1),mar=c(4,4.3,1,2.5))

plot(xy,xlim=c(-1.6,1.9),ylim=c(-1.9,1.3),bty="l",asp=1,xaxt="n",yaxt="n",type="n",ylab="",xlab="")
axis(side=1,at=seq(-2,2,0.5),cex.axis=0.9)
axis(side=2,at=seq(-2,2,0.5),las=2,cex.axis=0.9)
mtext(side=1, line=2.5,"dbRDA1 (37.38%)",cex=0.9)
mtext(side=2, line=3,"dbRDA2 (32.00%)",cex=0.9)

# Species names / scores
text(haps, rownames(haps),cex=0.6,col="darkgrey",xpd=T)

# Add location of each site.depth, but label by site only
# points(xy,pch=19,col=adjustcolor("white",alpha.f=0.6),cex=2.5)
text(xy, as.character(y_mat$Site),col=col.vec,cex=1)


# Add legend
legend(1,1.5, legend=xlabels,
       cex=0.7, bty="n",xpd=T,text.col=my_depth_symbols$color)
#######################################################################################









############# Supplementary Figure S1 #############
## Plot
xlabels=c("1m fringing","2m back reef","5m fore reef","10m fore reef","20m fore reef")
site.vec <- c(unique(raw_props$Site),"Pooled")

quartz(width=6,height=7)
par(mfrow=c(7,5),mar=c(1,0.5,0,0),oma=c(5,5,3,1))

ylims <- c(0.9,1,0.4,0.2,0.2,1)
bys <- c(0.3,0.3,0.1,0.1,0.1,0.3)

for(i in 1:7){
  for(j in 1:5){
    if(i!=7 & j!=5){
      y <- raw_props %>% filter(Species==cols$Species[i],Site==site.vec[j])
      y <- y[order(y$Depth.m),]
      x <- barplot(y$fit, yaxt="n", col=cols$Color[i], ylim=c(0,ylims[i]),border=NA)
      if(j==1) {axis(side=2,at=seq(0,ylims[i],bys[i]),las=1,cex.axis=1.2)}
      if(j!=1) {axis(side=2,at=seq(0,ylims[i],bys[i]),labels=NA)}
      axis(side=1,at=x,pos=-0.05*ylims[i],labels=NA)
      if(i==1){mtext(side=3,line=1,adj=0.5,eval(paste0("Site ",site.vec[j])))}
      # if(i==6){text(x,rep(-0.2,length(xlabels)),xlabels,srt=45,xpd=NA,adj=1,cex=1.2)}
    }
    if(i!=7 & j==5){
      y <- pred %>% filter(Species==cols$Species[i])
      y <- y[order(y$Depth.m),]
      x <- barplot(y$fit, yaxt="n", col=cols$Color[i], ylim=c(0,ylims[i]),border=NA)
      if(j==1) {axis(side=2,at=seq(0,ylims[i],bys[i]),las=1,cex.axis=1.2)}
      if(j!=1) {axis(side=2,at=seq(0,ylims[i],bys[i]),labels=NA)}
      segments(x,y$lwr,x,y$upr,col="black",lwd=1.5)
      axis(side=1,at=x,pos=-0.05*ylims[i],labels=NA)
      if(i==1){mtext(side=3,line=1,adj=0.5,eval(paste0("Sites ",site.vec[j])))}
    }
  if(i==7 & j!=5){
    y <- raw_props %>% filter(Site==site.vec[j])
    y <- y[order(y$Depth.m),]
    y_mat <- y %>% spread(Depth.m,fit)
    y_mat <- y_mat[match(cols$Species,y_mat$Species),] 
    y_mat[is.na(y_mat)] <- 0
    x <- barplot(as.matrix(y_mat[,-1:-2]),col=cols$Color,yaxt="n",xaxt="n",ylim=c(0,1.05),border=NA)
    axis(side=1,at=x,pos=-0.05,labels=NA)
    if(j==1) {axis(side=2,at=seq(0,1,0.3),las=1,cex.axis=1.2)}
    if(j!=1) {axis(side=2,at=seq(0,1,0.3),labels=NA)}
    text(x,rep(-0.2,length(xlabels)),xlabels,srt=45,xpd=NA,adj=1,cex=1.2)
  }
  if(i==7 & j==5){
    y_mat <- pred %>% select(Depth.m,Species,fit) %>% spread(Depth.m,fit)
    y_mat <- y_mat[match(cols$Species,y_mat$Species),] 
    y_mat[is.na(y_mat)] <- 0
    x <- barplot(as.matrix(y_mat[,-1]),col=cols$Color,yaxt="n",xaxt="n",ylim=c(0,1.05),border=NA)
    axis(side=1,at=x,pos=-0.05,labels=NA)
    text(x,rep(-0.2,length(xlabels)),xlabels,srt=45,xpd=NA,adj=1,cex=1.2)
  }
  }
}
mtext("Relative abundance (proportion of Pocillopora colonies)",side=2,outer=T,line=3,cex=1.2)

