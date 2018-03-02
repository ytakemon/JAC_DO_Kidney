######################################
# Analysis notes for Aging Kidney RNA/Protein project
# GAC 2-28-18
#####################################

library(tidyverse)
library(stringr)
library(magrittr)

######################################
# Load and format the data
######################################
load("./DO188b_kidney_201802_YT.Rdata")
#load("/ProcessedData/DO188b_kidney_201802_YT.Rdata")
ls()
# the data are in QTLViewer format
# we aren't doing mapping here so will clean up

###
#shorten names and remove genotype data
protein <- dataset.protein
rm(dataset.protein)
#
mrna <- dataset.mrna
rm(dataset.mrna)
#
rm(genoprobs, K, map, markers, genome.build)
#
ls()

####
# look at data structures

# the main data structure are lists: mrna and protein
names(protein)
names(mrna)

# covariates are in model matrix format
head(mrna$covar)
head(protein$covar)

# annotation keys are "id" which is Ens GeneID and Ens ProteinID
# the protein data includes GeneID that can be used to merge the structures
# both annotatiosn include gene "symbol"
head(mrna$annots)
head(protein$annots)

# in qtl2 style:  the column names are keys corresponding to "id"
colnames(protein$expr)[1:10]
colnames(mrna$expr)[1:10]

mrna$covar[mrna$covar[,3] == 1,][,3] <- "11"
mrna$covar[mrna$covar[,4] == 1,][,4] <- "12"
mrna$covar[mrna$covar[,5] == 1,][,5] <- "8"
mrna$covar[mrna$covar[,6] == 1,][,6] <- "9"

######################################
# Utility functions for pulling gene-level data
######################################
####
# get.gene()
# a function to extract data for one gene into a tibble
get.gene <- function(gene.symbol){

  # pull covariates and convert to factors
  # note Age is numeric and fAge is factor
  Mouse.ID <- rownames(mrna$covar)
  Sex <- factor(as.character(mrna$covar[,"SexM"]+1), levels=c("1","2"), labels=c("F","M"))
  Age <- mrna$covar[,"Age"]
  fAge <- factor(as.character(Age),levels=c("6","12","18"))
  G11 <- mrna$covar[,"GenerationG11"]
  G12 <- mrna$covar[,"GenerationG12"]
  G8 <- mrna$covar[,"GenerationG8"]
  G9 <- mrna$covar[,"GenerationG9"]

  #pull gene identifiers from mrna and protein annotations
  mrna.id <- mrna$annots[mrna$annots$symbol==gene.symbol,"id"]
  protein.id <- protein$annots[protein$annots$symbol==gene.symbol,"id"]

  # strat building the tiblle to hold data
  my.data <- data_frame(Mouse.ID=Mouse.ID,
                        Sex=Sex,
                        Age=Age,
                        fAge=fAge,
                        G11=G11,
                        G12=G12,
                        G8=G8,
                        G9=G9)
  # add RNA and Protein columns to the tibble - if it exists
  if(length(mrna.id)>0){
    my.data <- mutate(my.data, RNA=mrna$expr[,mrna.id])
  }
  if(length(protein.id)>0){
    my.data <- mutate(my.data, Protein=protein$expr[,protein.id])
  }

  # return value is a tibble with covariates, RNA and Protein data
  my.data
}
# #test
# get.gene("Gnai3")
# get.gene("Cdkn2a")  # a case with no Protein data
# get.gene("Cow")      # a non-existing gene

###
# how to make a boxplot of gene data
quartz()
get.gene("Cdkn2a") %>%
  ggplot(aes(x=Sex:fAge, y=RNA, color=Sex)) +
      geom_boxplot() +
      geom_point(position = position_jitter(width=0.25, height=0)) +
      ggtitle("Cdkn2a")

###
# how to make a scatter plot of gene data
quartz()
get.gene("Gnai3") %>%
  ggplot(aes(x=RNA, y=Protein, color=fAge)) +
    geom_point() +
    geom_smooth(method="lm") +
    facet_wrap(~ Sex)

###
# how to fit a linear model to gene data
anova(lm(Protein ~ Sex+fAge+Sex:fAge+RNA+Sex:RNA+fAge:RNA,
         data=get.gene("Gnai3")))

######################################
# Create a data summary with GeneID as key
#
# recreate Yuka's figure 1
######################################

###
# create a tibble with covariates
covs <- get.gene("Cow")

###
# function to get slopes
get.slope <- function(df){
  coefficients(lm(GEX ~ Sex+Age, data=df))["Age"]
}

###
# first summarize the mRNA data

# regression slopes
mrna.slope <- as_tibble(mrna$expr) %>%
  mutate(Sex=covs$Sex, Age=covs$Age) %>%
  gather(starts_with("ENS"),key="GeneID",value="GEX") %>%
  group_by(GeneID) %>%
  nest() %>%
  mutate(rna.slope=map_dbl(data, get.slope)) %>%
  select(GeneID, rna.slope)

# mrna mean expression by groups
mrna.means <- as_tibble(mrna$expr) %>%
  mutate(Sex=covs$Sex, fAge=covs$fAge) %>%
  gather(starts_with("ENS"),key="GeneID",value="RNA") %>%
  group_by(Sex, fAge, GeneID) %>%
  summarize(MeanRNA=mean(RNA)) %>%
  unite(Group, Sex, fAge, sep="_") %>%
  spread(key=Group,value=MeanRNA) %>%
  rename(rna_F_6=F_6,rna_F_12=F_12,rna_F_18=F_18,
         rna_M_6=M_6,rna_M_12=M_12,rna_M_18=M_18) %>%
  left_join(mrna.slope)
###
# repeat for protein data

# regression slopes
protein.slope <- as_tibble(protein$expr) %>%
  mutate(Sex=covs$Sex, Age=covs$Age) %>%
  gather(starts_with("ENS"),key="ProteinID",value="GEX") %>%
  group_by(ProteinID) %>%
  nest() %>%
  mutate(protein.slope=map_dbl(data, get.slope)) %>%
  select(ProteinID, protein.slope)

# protein mean expression by groups
protein.means <- as_tibble(protein$expr) %>%
  mutate(Sex=covs$Sex, fAge=covs$fAge) %>%
  gather(starts_with("ENS"),key="ProteinID",value="GEX") %>%
  group_by(Sex, fAge, ProteinID) %>%
  summarize(MeanProtein=mean(GEX, na.rm=TRUE)) %>%
  unite(Group, Sex, fAge, sep="_") %>%
  spread(key=Group,value=MeanProtein) %>%
  rename(protein_F_6=F_6, protein_F_12=F_12, protein_F_18=F_18,
         protein_M_6=M_6, protein_M_12=M_12, protein_M_18=M_18) %>%
  left_join(protein.slope)


##########
# Figure out differences between Gary's and my calculation
#########

# Adding generation
get.slope.test <- function(df){
  coefficients(lm(GEX ~ Sex+Age+G8+G9+G11+G12, data=df))["Age"]
}

mrna.slope.test <- as_tibble(mrna$expr) %>%
  mutate(Sex=covs$Sex, Age=covs$Age, G8=covs$G8, G9=covs$G9, G11=covs$G11, G12=covs$G12) %>%
  gather(starts_with("ENS"),key="GeneID",value="GEX") %>%
  group_by(GeneID) %>%
  nest() %>%
  mutate(rna.slope=map_dbl(data, get.slope.test)) %>%
  select(GeneID, rna.slope)

protein.slope.test <- as_tibble(protein$expr) %>%
  mutate(Sex=covs$Sex, Age=covs$Age, G8=covs$G8, G9=covs$G9, G11=covs$G11, G12=covs$G12) %>%
  gather(starts_with("ENS"),key="ProteinID",value="GEX") %>%
  group_by(ProteinID) %>%
  nest() %>%
  mutate(protein.slope=map_dbl(data, get.slope.test)) %>%
  select(ProteinID, protein.slope)

protein.id <- as_tibble(protein$annots) %>%
  select(ProteinID=id, GeneID=gene_id)
protein.slope.test<- left_join(protein.id, protein.slope.test)

slope.test <- inner_join(mrna.slope.test, protein.slope.test, by="GeneID")

ggplot(slope.test, aes(x=rna.slope, y=protein.slope)) +
  geom_point()
###########
##########


#add ensembl GeneID
protein.id <- as_tibble(protein$annots) %>%
  select(ProteinID=id, GeneID=gene_id)
#
protein.means <- left_join(protein.id, protein.means)

####
# join the mRNA and protein data for shared GeneIDs
expr.means <- inner_join(mrna.means, protein.means, by="GeneID")
#
dim(expr.means)
# there are 6667 gene with both mRNA and protein data

####
#clean up
rm(mrna.slope, mrna.means, protein.slope, protein.means, protein.id)

####
quartz()
ggplot(expr.means, aes(x=rna.slope, y=protein.slope)) +
  geom_point()







######################################
# Slc9a3
######################################

####
# look at the full anova model for rgression of Protein on RNA
anova(lm(Protein ~ Sex+fAge+Sex:fAge+RNA+Sex:RNA+fAge:RNA, data=get.gene("Slc9a3")))
# the only significant terms are Sex and fAge

####
# look at the regression slope or Protein on RNA in reduced model
coefficients(lm(Protein ~ Sex+fAge+RNA, data=get.gene("Slc9a3")))
# slope is 0.07 (N.S.)

####
# build up Slc9a3 data for plotting
Slc9a3 <- get.gene("Slc9a3")

#add fitted regression line for plotting
Slc9a3 <- mutate(Slc9a3,
                 Fitted = fitted(lm(Protein ~ Sex+fAge+Sex:fAge+RNA,
                                    data=Slc9a3)))

####
# summarize means and sd for plotting
Slc9a3.summ <- Slc9a3 %>%
  group_by(Sex, fAge) %>%
  summarize(mRNA=mean(RNA), sdRNA=sd(RNA),
            mProtein=mean(Protein), sdProtein=sd(Protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

####
# the default plot with seperate regression by groups
quartz()
get.gene("Slc9a3") %>%
  ggplot(aes(x=RNA, y=Protein, color=fAge)) +
  geom_point() +
  geom_smooth(method="lm") +
  facet_wrap(~ Sex)

####
# the plot with common slope from regression fit
quartz()
ggplot(Slc9a3.summ, aes(x=mRNA, y=mProtein, color=fAge, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=Slc9a3, aes(x=RNA, y=Protein)) +
  geom_line(data=Slc9a3, aes(x = RNA, y = Fitted)) +
  ggtitle("Slc9a3")

#####
# INCOMPLETE
# write a function to make this pplot
my.plot1 <- function(gene.name){
  EnsemblID <- annot.rna[annot.rna$Gene==gene.name,"EnsemblID"]
  if(length(EnsemblID)>0){
  p <- as_tibble(expr.z) %>%
   select_(EnsemblID) %>%
    mutate(Mouse.ID=Mouse.ID, Group=Sex:Diet) %>%
    ggplot(aes_string(x="Group", y=EnsemblID, color="Group")) +
      geom_boxplot() +
      geom_point(position = position_jitter(width=0.25, height=0)) +
      ggtitle(gene.name)
  }
  else{ p <- gene.name}
  p
}
