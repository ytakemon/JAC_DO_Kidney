setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
library(tidyverse)
library(stringr)

######################################
# Load and format the data
######################################
load("./RNAseq_data/DO188b_kidney_201802_YT.Rdata")
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

  #pull gene identifiers from mrna and protein annotations
  mrna.id <- mrna$annots[mrna$annots$symbol==gene.symbol,"id"]
  protein.id <- protein$annots[protein$annots$symbol==gene.symbol,"id"]

  # strat building the tiblle to hold data
  my.data <- data_frame(Mouse.ID=Mouse.ID,
                        Sex=Sex,
                        Age=Age,
                        fAge=fAge)
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

## INCOMPLETE ###
######################################

###
# create a tibble with covariates
covs <- get.gene("Cow")

###
# first summarize the mRNA data

# regression slopes
mrna.slope <- as_tibble(mrna$expr) %>%
  mutate(Sex=covs$Sex, Age=covs$Age) %>%
  gather(starts_with("ENS"),key="GeneID",value="RNA") %>%
  group_by(GeneID) %>%
  summarize(rna.slope = coefficients(lm(RNA ~ Sex+Age, data=genedata))["Age"])

# mrna mean expression by groups
mrna.means <- as_tibble(mrna$expr) %>%
  mutate(Sex=covs$Sex, fAge=covs$fAge) %>%
  gather(starts_with("ENS"),key="GeneID",value="RNA") %>%
  group_by(Sex, fAge, GeneID) %>%
  summarize(MeanRNA=mean(RNA)) %>%
  unite(Group, Sex, fAge, sep="_") %>%
  spread(key=Group,value=MeanRNA)


%>%
  left_join(as_tibble(annot.rna) %>%
                      select(EnsemblID, Gene)) %>%
  select(EnsemblID, Gene, everything())

#look
expr.t


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
