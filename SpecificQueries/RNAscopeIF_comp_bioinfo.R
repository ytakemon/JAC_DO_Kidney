# Identify animals that have the average given expression values for the given group.
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggsci)
library(stringr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")
load("./shiny_annotation.RData")

# Selection quality:
Gene_name <- "Lrp2"
Data <- "protein"
Age <- c(6, 18)

samples <- annot.samples
# Identify gene name
other.ids <- function(gene.name, level) {
  if (level == "mRNA") {
    sel <- which(mRNA.list$symbol == gene.name)[1]
    if (!is.na(sel)) return(mRNA.list[sel,]) else return(c(NA,NA,NA))
  }
  if (level == "protein") {
    sel <- which(protein.list$symbol == gene.name)[1]
    if (!is.na(sel)) return(protein.list[sel,]) else return(c(NA,NA,NA))
  }
}

gene <- other.ids(Gene_name, Data)
# Get expression data
samples$RNAexpr <- expr.mrna[, gene$id]
samples$Protexpr <- expr.protein[, gene$protein_id]

# Get RNAscope data
df <- read.csv("./Phenotype/phenotypes/Lrp2RNAscope_cleaned.csv")
df$ratio <- df$RNAscope_Lrp2_count / df$IF_LRP2_intensity

# Take log of genometic mean for each sample
animal <- as.data.frame(unique(df$Sample))
animal$Age <- NA
# Get age
for (i in 1:nrow(animal)){
  animal$Age[i] <- df[df$Sample == animal[,1][i],]$Age[1]
}
animal$log_mean_ratio <- NA
animal$log_mean_rna <- NA
animal$log_mean_protein <- NA
# Get geometic mean for each sample
for (i in 1:nrow(animal)){
  sub <- df[df$Sample == animal[,1][i],]
  animal$log_mean_ratio[i] <- log(mean(sub$ratio))
  animal$log_mean_rna[i] <- log(mean(sub$RNAscope_Lrp2_count))
  animal$log_mean_protein[i] <- log(mean(sub$IF_LRP2_intensity))
}
animal <- arrange(animal, Age)
names(animal)[1] <- "Sample"
animal$Age <- as.factor(as.character(animal$Age))
animal$Age <- factor(animal$Age,levels(animal$Age)[c(2,1)])
animal <- arrange(animal, Sample)
rownames(animal) <- animal$Sample

# Subset RNAseq/ Proteomics data
samples <- samples[rownames(samples) %in% rownames(animal),]
samples <- arrange(samples, Mouse.ID)
samples$log_mean_ratio <- animal$log_mean_ratio
samples$log_mean_rna <- animal$log_mean_rna
samples$log_mean_protein <- animal$log_mean_protein
samples$ratio <- samples$RNAexpr / samples$Protexpr

ggplot(samples, aes(x = ratio, y = log_mean_ratio, colour = as.factor(as.character(Age)))) +
  geom_point() +
  geom_smooth( method = lm, se = FALSE) +
  labs( title = "") +
  scale_color_aaas()

ggplot(samples, aes(x = RNAexpr, y = Protexpr, colour = as.factor(as.character(Age)))) +
  geom_point() +
  geom_smooth( method = lm, se = FALSE) +
  labs( title = "") +
  scale_color_aaas()

ggplot(samples, aes(x = log_mean_rna, y = log_mean_protein, colour = as.factor(as.character(Age)))) +
  geom_point() +
  geom_smooth( method = lm, se = FALSE) +
  labs( title = "") +
  scale_color_aaas()

ggplot(samples, aes(x = RNAexpr, y = log_mean_rna, colour = as.factor(as.character(Age)))) +
  geom_point() +
  geom_smooth( method = lm, se = FALSE) +
  labs( title = "", x = "RNA-seq Lrp2", y = "RNAscopr Lrp2") +
  scale_color_aaas()

ggplot(samples, aes(x = Protexpr, y = log_mean_protein, colour = as.factor(as.character(Age)))) +
  geom_point() +
  geom_smooth( method = lm, se = FALSE) +
  labs( title = "", x = "Proteomics Lrp2", y = "IF Lrp2") +
  scale_color_aaas()
