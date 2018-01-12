# Identify animals that have the average given expression values for the given group.
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggsci)
library(stringr)
library(gridExtra)
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
samples$RNAraw <- raw.mrna[, gene$id]
samples$Protraw <- raw.protein[, gene$protein_id]

# Get RNAscope data v1 -------------------------------------------------------
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
  animal$log_mean_ratio[i] <- log(mean(sub$ratio, na.rm = TRUE))
  animal$log_mean_rna[i] <- log(mean(sub$RNAscope_Lrp2_count, na.rm = TRUE))
  animal$log_mean_protein[i] <- log(mean(sub$IF_LRP2_intensity, na.rm = TRUE))
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
samples$Age <- as.factor(as.character(samples$Age))

# RNAscope v2-----------------------------------------------------------------
df <- read.csv("./Phenotype/phenotypes/Lrp2RNAscope_trial2.csv")
df <- df[,c(1:6,11,12)]

# calculate
df$ProteinRNAratio <- df$Protein / df$CountRNA
df$ProteinRNAratio[df$ProteinRNAratio %in% c(Inf, NA)] <- NA
df$ProteinRNAratioPerTubArea <- df$ProteinRNAratio / df$Tubule.Area
df$ProteinRNAratioPerTubArea[df$ProteinRNAratioPerTubArea %in% c(Inf, NA)] <- NA
df$ProteinRNAratioPerNuclei <- df$ProteinRNAratio / df$X.Nuclei
df$ProteinRNAratioPerNuclei[df$ProteinRNAratioPerNuclei %in% c(Inf, NA)] <- NA
df$CountRNA[df$CountRNA %in% c(Inf, NA)] <- NA
df$Protein[df$Protein %in% c(Inf, NA)] <- NA

# Make new df and get age
animal <- as.data.frame(as.character(unique(df$Mouse)))
animal$Age <- NA
for (i in 1:nrow(animal)){
  animal$Age[i] <- as.character(df[df$Mouse == animal[,1][i],]$Group[1])
}
animal <- arrange(animal, Age)
names(animal)[1] <- "Sample"

# Calc log - geometric mean
animal$log_mean_ratio <- NA
animal$log_mean_ratioPerTub <- NA
animal$log_mean_ratioPerNuc <- NA
animal$log_mean_RNA <- NA
animal$log_mean_protein <- NA

for (i in 1:nrow(animal)){
  sub <- df[df$Mouse == animal[,1][i],]
  animal$log_mean_ratio[i] <- log(mean(sub$ProteinRNAratio, na.rm = TRUE))
  animal$log_mean_ratioPerTub[i] <- log(mean(sub$ProteinRNAratioPerTubArea, na.rm = TRUE))
  animal$log_mean_ratioPerNuc[i] <- log(mean(sub$ProteinRNAratioPerNuclei, na.rm = TRUE))
  animal$log_mean_RNA[i] <- log(mean(sub$CountRNA, na.rm = TRUE))
  animal$log_mean_protein[i] <- log(mean(sub$Protein, na.rm = TRUE))
}

# make factor and reorder levels
animal$Age <- as.factor(as.character(animal$Age))
animal$Age <- factor(animal$Age, levels(animal$Age)[c(2,1)])
animal$Sample <- paste0("DO-",animal$Sample)
animal <- arrange(animal, Sample)

#Fix typos......again.
animal$Sample[5] <- "DO-1263"
animal$Sample[7] <- "DO-0937"
animal$Sample[8] <- "DO-0950"
animal$Sample[9] <- "DO-0952"

animal <- arrange(animal, Sample)
samples <- arrange(samples, Mouse.ID)

# Plot -------------------------------------------------------------------------
plotDir <- "/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Plot"

RNArank <- ggplot(samples, aes(x = RNAexpr, y = log_mean_rna, colour = Age)) +
  geom_point() +
  geom_smooth( method = lm, se = FALSE) +
  labs( title = "RANK RNA", x = "RNA-seq Lrp2 (rank norm)", y = "RNAscopr Lrp2 (log norm)") +
  scale_color_aaas()

Protrank <- ggplot(samples, aes(x = Protexpr, y = log_mean_protein, colour = Age)) +
  geom_point() +
  geom_smooth( method = lm, se = FALSE) +
  labs( title = "RANK PROT", x = "Proteomics Lrp2 (rank norm)", y = "IF Lrp2 (log norm)") +
  scale_color_aaas()

RNAraw <- ggplot(samples, aes(x = log(RNAraw), y = log_mean_rna, colour = Age)) +
  geom_point() +
  geom_smooth( method = lm, se = FALSE) +
  labs( title = "LOG RAW RNA", x = "RNA-seq Lrp2 (log raw)", y = "RNAscopr Lrp2 (log norm)") +
  scale_color_aaas()

Protraw <- ggplot(samples, aes(x = log(Protraw), y = log_mean_protein, colour = Age)) +
  geom_point() +
  geom_smooth( method = lm, se = FALSE) +
  labs( title = "LOG RAW PROT", x = "Proteomics Lrp2 (log raw)", y = "IF Lrp2 (log norm)") +
  scale_color_aaas()

pdf(paste0(plotDir,"/RankCompare_Lpr2_RNAv1.pdf"), height = 6, width = 12)
grid.arrange(RNArank, Protrank, ncol = 2)
dev.off()

pdf(paste0(plotDir,"/RawCompare_Lpr2_RNAv1.pdf"), height = 6, width = 12)
grid.arrange(RNAraw, Protraw, ncol = 2)
dev.off()

# Plot V2 ----------------------------------------------------------------------

# add v2 data
samples <- samples[samples$Mouse.ID %in% animal$Sample,]
samples$log_mean_rna_v2 <- animal$log_mean_RNA
samples$log_mean_prot_v2 <- animal$log_mean_protein

RNAraw <- ggplot(samples, aes(x = log(RNAraw), y = log_mean_rna_v2, colour = Age)) +
  geom_point() +
  geom_smooth( method = lm, se = FALSE) +
  labs( title = "V2 LOG RAW RNA", x = "RNA-seq Lrp2 (log raw)", y = "RNAscopr Lrp2 (log norm)") +
  scale_color_aaas()

Protraw <- ggplot(samples, aes(x = log(Protraw), y = log_mean_prot_v2, colour = Age)) +
  geom_point() +
  geom_smooth( method = lm, se = FALSE) +
  labs( title = "V2 LOG RAW PROT", x = "Proteomics Lrp2 (log raw)", y = "IF Lrp2 (log norm)") +
  scale_color_aaas()

pdf(paste0(plotDir,"/RawCompare_Lpr2_RNAv2.pdf"), height = 6, width = 12)
grid.arrange(RNAraw, Protraw, ncol = 2)
dev.off()
