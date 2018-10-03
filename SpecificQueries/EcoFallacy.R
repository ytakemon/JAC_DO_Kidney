######################################
# Analysis notes for Aging Kidney RNA/Protein project
# GAC 2-28-18
#####################################

library(tidyverse)
library(stringr)
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
df <- get.gene("Slc9a3") %>%
  mutate(
    Fitted = fitted(lm(Protein ~ Sex+fAge+Sex:fAge+RNA, data=.)))

df_summary <- df %>%
  group_by(Sex, fAge) %>%
  summarize(mRNA=mean(RNA), sdRNA=sd(RNA),
    mProtein=mean(Protein), sdProtein=sd(Protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

pdf("./Plot/Slc9a3cor_byAge.pdf", width = 8, height = 7)
ggplot(df_summary, aes(x=mRNA, y=mProtein, color=fAge, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df, aes(x=RNA, y=Protein)) +
  geom_line(data=df, aes(x = RNA, y = Fitted)) +
  labs( title = "Slc9a3 mRNA v. Protein Expression by Age",
        subtitle = paste0("Female: 93 (6mo=33, 12mo=31, 18mo=29) \nMale: 95(6mo=30, 12mo=31, 18mo=34)"),
        y = "SLC9A3 expression (Rank normalized protein)",
        x = "Slc9a3 expression (Rank normalized mRNA)",
        color = "Age") +
  theme_bw() +
  scale_colour_aaas()
dev.off()

######################################
# Gprc5c
######################################

####
# look at the full anova model for rgression of Protein on RNA
anova(lm(Protein ~ Sex+fAge+Sex:fAge+RNA+Sex:RNA+fAge:RNA, data=get.gene("Gprc5c")))
# the only significant terms are Sex and fAge and RNA

####
# look at the regression slope or Protein on RNA in reduced model
coefficients(lm(Protein ~ Sex+fAge+RNA, data=get.gene("Gprc5c")))
# slope is NS

####
# build up Gprc5c data for plotting
df <- get.gene("Gprc5c") %>%
  mutate(
    Fitted = fitted(lm(Protein ~ Sex+fAge+Sex:fAge+RNA, data=.)))

df_summary <- df %>%
  group_by(Sex, fAge) %>%
  summarize(mRNA=mean(RNA), sdRNA=sd(RNA),
    mProtein=mean(Protein), sdProtein=sd(Protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

pdf("./Plot/Gprc5ccor_byAge.pdf", width = 8, height = 7)
ggplot(df_summary, aes(x=mRNA, y=mProtein, color=fAge, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df, aes(x=RNA, y=Protein)) +
  geom_line(data=df, aes(x = RNA, y = Fitted)) +
  labs( title = "Gprc5c mRNA v. Protein Expression by Age",
        subtitle = paste0("Female: 93 (6mo=33, 12mo=31, 18mo=29) \nMale: 95(6mo=30, 12mo=31, 18mo=34)"),
        y = "GPRC5C expression (Rank normalized protein)",
        x = "Gprc5c expression (Rank normalized mRNA)",
        color = "Age") +
  theme_bw() +
  scale_colour_aaas()
dev.off()

######################################
# Slc5a12
######################################

####
# look at the full anova model for rgression of Protein on RNA
anova(lm(Protein ~ Sex+fAge+Sex:fAge+RNA+Sex:RNA+fAge:RNA, data=get.gene("Slc5a12")))
# the only significant terms are fAge and RNA

####
# look at the regression slope or Protein on RNA in reduced model
coefficients(lm(Protein ~ Sex+fAge+RNA, data=get.gene("Slc5a12")))
# slope is NS

####
# build up Slc5a12 data for plotting
df <- get.gene("Slc5a12") %>%
  mutate(
    Fitted = fitted(lm(Protein ~ Sex+fAge+Sex:fAge+RNA, data=.)))

df_summary <- df %>%
  group_by(Sex, fAge) %>%
  summarize(mRNA=mean(RNA), sdRNA=sd(RNA),
    mProtein=mean(Protein), sdProtein=sd(Protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

pdf("./Plot/Slc5a12cor_byAge.pdf", width = 8, height = 7)
ggplot(df_summary, aes(x=mRNA, y=mProtein, color=fAge, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df, aes(x=RNA, y=Protein)) +
  geom_line(data=df, aes(x = RNA, y = Fitted)) +
  labs( title = "Slc5a12 mRNA v. Protein Expression by Age",
        subtitle = paste0("Female: 93 (6mo=33, 12mo=31, 18mo=29) \nMale: 95(6mo=30, 12mo=31, 18mo=34)"),
        y = "SLC5A12 expression (Rank normalized protein)",
        x = "Slc5a12 expression (Rank normalized mRNA)",
        color = "Age") +
  theme_bw() +
  scale_colour_aaas()
dev.off()

######################################
# Erh
######################################

####
# look at the full anova model for rgression of Protein on RNA
anova(lm(Protein ~ Sex+fAge+Sex:fAge+RNA+Sex:RNA+fAge:RNA, data=get.gene("Erh")))
# the only significant terms are fAge, Sex:fAge, and fAge:RNA

####
# look at the regression slope or Protein on RNA in reduced model
coefficients(lm(Protein ~ fAge + Sex:fAge+fAge:RNA, data=get.gene("Erh")))
# slope is -0.97563818 at fAge18

####
# build up Erh data for plotting
df <- get.gene("Erh")
df <- df[complete.cases(df$Protein),]%>%
  mutate(
    Fitted = fitted(lm(Protein ~ Sex+fAge+Sex:fAge+RNA, data=.)))

df_summary <- df %>%
  group_by(Sex, fAge) %>%
  summarize(mRNA=mean(RNA), sdRNA=sd(RNA),
    mProtein=mean(Protein), sdProtein=sd(Protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

pdf("./Plot/Erhcor_byAge.pdf", width = 8, height = 7)
ggplot(df_summary, aes(x=mRNA, y=mProtein, color=fAge, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df, aes(x=RNA, y=Protein)) +
  geom_line(data=df, aes(x = RNA, y = Fitted)) +
  labs( title = "Erh mRNA v. Protein Expression by Age",
        subtitle = paste0("Female: 93 (6mo=33, 12mo=31, 18mo=29) \nMale: 95(6mo=30, 12mo=31, 18mo=34)"),
        y = "ERH expression (Rank normalized protein)",
        x = "Erh expression (Rank normalized mRNA)",
        color = "Age") +
  theme_bw() +
  scale_colour_aaas()
dev.off()


######################################
# Flot1
######################################

####
# look at the full anova model for rgression of Protein on RNA
anova(lm(Protein ~ Sex+fAge+Sex:fAge+RNA+Sex:RNA+fAge:RNA, data=get.gene("Flot1")))
# the only significant terms are fAge and RNA
####
# look at the regression slope or Protein on RNA in reduced model
coefficients(lm(Protein ~ fAge + RNA, data=get.gene("Flot1")))
# slope is -1.6657549 at fAge18

####
# build up Flot1 data for plotting
df <- get.gene("Flot1") %>%
  mutate(
    Fitted = fitted(lm(Protein ~ Sex+fAge+Sex:fAge+RNA, data=.)))

df_summary <- df %>%
  group_by(Sex, fAge) %>%
  summarize(mRNA=mean(RNA), sdRNA=sd(RNA),
    mProtein=mean(Protein), sdProtein=sd(Protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

pdf("./Plot/Flot1cor_byAge.pdf", width = 8, height = 7)
ggplot(df_summary, aes(x=mRNA, y=mProtein, color=fAge, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df, aes(x=RNA, y=Protein)) +
  geom_line(data=df, aes(x = RNA, y = Fitted)) +
  labs( title = "Flot1 mRNA v. Protein Expression by Age",
        subtitle = paste0("Female: 93 (6mo=33, 12mo=31, 18mo=29) \nMale: 95(6mo=30, 12mo=31, 18mo=34)"),
        y = "FLOT1 expression (Rank normalized protein)",
        x = "Flot1 expression (Rank normalized mRNA)",
        color = "Age") +
  theme_bw() +
  scale_colour_aaas()
dev.off()


######################################
# Rnf121
######################################

####
# look at the full anova model for rgression of Protein on RNA
anova(lm(Protein ~ Sex+fAge+Sex:fAge+RNA+Sex:RNA+fAge:RNA, data=get.gene("Rnf121")))
# the only significant terms is fAge
####
# look at the regression slope or Protein on RNA in reduced model
coefficients(lm(Protein ~ fAge, data=get.gene("Rnf121")))
# slope is -1.1681476 at fAge18

####
# build up Rnf121 data for plotting
df <- get.gene("Rnf121")
df <- df[complete.cases(df$Protein),]%>%
  mutate(
    Fitted = fitted(lm(Protein ~ Sex+fAge+Sex:fAge+RNA, data=.)))

df_summary <- df %>%
  group_by(Sex, fAge) %>%
  summarize(mRNA=mean(RNA), sdRNA=sd(RNA),
    mProtein=mean(Protein), sdProtein=sd(Protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

pdf("./Plot/Rnf121cor_byAge.pdf", width = 8, height = 7)
ggplot(df_summary, aes(x=mRNA, y=mProtein, color=fAge, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df, aes(x=RNA, y=Protein)) +
  geom_line(data=df, aes(x = RNA, y = Fitted)) +
  labs( title = "Rnf121 mRNA v. Protein Expression by Age",
        subtitle = paste0("Female: 93 (6mo=33, 12mo=31, 18mo=29) \nMale: 95(6mo=30, 12mo=31, 18mo=34)"),
        y = "RNF121 expression (Rank normalized protein)",
        x = "Rnf121 expression (Rank normalized mRNA)",
        color = "Age") +
  theme_bw() +
  scale_colour_aaas()
dev.off()

####
# Combine Slc5a12 and Flot1 into one figure
# Figure 3A is Slc5a12 (top)
# Fiugure 3B is Flot1 (bottom)

# Slc5a12
df <- get.gene("Slc5a12") %>%
  mutate(
    Fitted = fitted(lm(Protein ~ Sex+fAge+Sex:fAge+RNA, data=.)))

df_summary <- df %>%
  group_by(Sex, fAge) %>%
  summarize(mRNA=mean(RNA), sdRNA=sd(RNA),
    mProtein=mean(Protein), sdProtein=sd(Protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

Slc5a12 <- ggplot(df_summary, aes(x=mRNA, y=mProtein, color=fAge, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df, aes(x=RNA, y=Protein)) +
  geom_line(data=df, aes(x = RNA, y = Fitted)) +
  labs( title = "Slc5a12 mRNA v. Protein Expression by Age",
        subtitle = paste0("Female: 93 (6mo=33, 12mo=31, 18mo=29) \nMale: 95(6mo=30, 12mo=31, 18mo=34)"),
        y = "SLC5A12 expression (Rank normalized protein)",
        x = "Slc5a12 expression (Rank normalized mRNA)",
        color = "Age") +
  theme_bw() +
  scale_colour_aaas()

# Flot1
df <- get.gene("Flot1") %>%
  mutate(
    Fitted = fitted(lm(Protein ~ Sex+fAge+Sex:fAge+RNA, data=.)))

df_summary <- df %>%
  group_by(Sex, fAge) %>%
  summarize(mRNA=mean(RNA), sdRNA=sd(RNA),
    mProtein=mean(Protein), sdProtein=sd(Protein), N=n() ) %>%
  mutate(seRNA = sdRNA/sqrt(N), seProtein = sdProtein/sqrt(N))

Flot1 <- ggplot(df_summary, aes(x=mRNA, y=mProtein, color=fAge, shape=Sex), size=4) +
  geom_point(size=4) +
  geom_errorbar(aes(ymin=mProtein-seProtein, ymax=mProtein+seProtein),size=1) +
  geom_errorbarh(aes(xmin=mRNA-seRNA, xmax=mRNA+seRNA),size=1) +
  geom_point(data=df, aes(x=RNA, y=Protein)) +
  geom_line(data=df, aes(x = RNA, y = Fitted)) +
  labs( title = "Flot1 mRNA v. Protein Expression by Age",
        subtitle = paste0("Female: 93 (6mo=33, 12mo=31, 18mo=29) \nMale: 95(6mo=30, 12mo=31, 18mo=34)"),
        y = "FLOT1 expression (Rank normalized protein)",
        x = "Flot1 expression (Rank normalized mRNA)",
        color = "Age") +
  theme_bw() +
  scale_colour_aaas()


library(gridExtra)
pdf("./Plot/Figure3_EcoFallacy.pdf", width = 5, height = 10)
grid.arrange(Slc5a12, Flot1, ncol =1)
dev.off()
