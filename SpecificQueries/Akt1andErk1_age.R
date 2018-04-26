######################################
# Analysis notes for Aging Kidney RNA/Protein project
# GAC 2-28-18
#####################################

library(tidyverse)
library(stringr)
library(magrittr)
library(ggsci)

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

######################################
# Akt1
######################################

####
# look at the full anova model for rgression of Protein on RNA
anova(lm(Protein ~ Sex+fAge+Sex:fAge+RNA+Sex:RNA+fAge:RNA, data=get.gene("Akt1")))
# the only significant Sex, fAge, RNA

####
# look at the regression slope or Protein on RNA in reduced model
coefficients(lm(Protein ~ Sex+fAge+RNA, data=get.gene("Akt1")))
# slope is 0.9 at fAge18

####
# build up Akt1 data for plotting
df <- get.gene("Akt1") %>%
  mutate(
    Fitted = fitted(lm(Protein ~ Sex+fAge+Sex:fAge+RNA, data=.)))

df_summary <- df %>%
  group_by(Sex, fAge) %>%
  summarize(mRNA=mean(RNA), sdRNA=sd(RNA),
    mProtein=mean(Protein), sdProtein=sd(Protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

pdf("./Plot/Akt1cor_byAge.pdf", width = 8, height = 7)
ggplot(df_summary, aes(x=mRNA, y=mProtein, color=fAge, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df, aes(x=RNA, y=Protein)) +
  geom_line(data=df, aes(x = RNA, y = Fitted)) +
  labs( title = "Akt1 mRNA v. Protein Expression by Age",
        subtitle = paste0("Female: 93 (6mo=33, 12mo=31, 18mo=29) \nMale: 95(6mo=30, 12mo=31, 18mo=34)"),
        y = "AKT1 expression (Rank normalized protein)",
        x = "Akt1 expression (Rank normalized mRNA)",
        color = "Age") +
  theme_bw() +
  scale_colour_aaas()
dev.off()

######################################
# Erk1
######################################

####
# look at the full anova model for rgression of Protein on RNA
anova(lm(Protein ~ Sex+fAge+Sex:fAge+RNA+Sex:RNA+fAge:RNA, data=get.gene("Mapk3")))
# the only significant Sex, fAge, RNA

####
# look at the regression slope or Protein on RNA in reduced model
coefficients(lm(Protein ~ Sex+fAge+RNA, data=get.gene("Mapk3")))
# slope is 1.04 at fAge18

####
# build up Mapk3 data for plotting
df <- get.gene("Mapk3") %>%
  mutate(
    Fitted = fitted(lm(Protein ~ Sex+fAge+Sex:fAge+RNA, data=.)))

df_summary <- df %>%
  group_by(Sex, fAge) %>%
  summarize(mRNA=mean(RNA), sdRNA=sd(RNA),
    mProtein=mean(Protein), sdProtein=sd(Protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

pdf("./Plot/Erk1cor_byAge.pdf", width = 8, height = 7)
ggplot(df_summary, aes(x=mRNA, y=mProtein, color=fAge, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df, aes(x=RNA, y=Protein)) +
  geom_line(data=df, aes(x = RNA, y = Fitted)) +
  labs( title = "Erk1 mRNA v. Protein Expression by Age",
        subtitle = paste0("Female: 93 (6mo=33, 12mo=31, 18mo=29) \nMale: 95(6mo=30, 12mo=31, 18mo=34)"),
        y = "ERK1 expression (Rank normalized protein)",
        x = "Erk1 expression (Rank normalized mRNA)",
        color = "Age") +
  theme_bw() +
  scale_colour_aaas()
dev.off()
