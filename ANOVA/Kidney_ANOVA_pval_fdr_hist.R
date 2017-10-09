# ANOVA for both RNA and protein:
# 1. Age (adjusted for sex and generation)
# 2. Sex (adjusted for age and generation)
# 3. Age:Sex interaction with age as a variate (adjust for generation).
# (Compute p-values and FDR using BH method)

# References:
# https://github.com/simecek/TheAgingProteome/blob/master/analysis/pval-histograms.Rmd

# Usage on Cadillac
# cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts
# qsub -v I=kidney_anova_fdr_output.csv,script=Kidney_ANOVA_pval_hist Rsubmit_args.sh

# R/3.3.2
library(ggplot2)
library(gridExtra)

# Load data
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
args <- commandArgs(trailingOnly = TRUE)
output <- list.files(pattern = paste0("^",args[[1]]), recursive = TRUE)
output <- read.csv(output[[1]])

# Plot pvalues
# 1.1.1 mRNA - Age (adjusted for sex and generation)
pval_mRNA_age <- ggplot(output, aes(x=p.mRNA_Age.Sex)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="mRNA/Age",
         x="p-value")
# 1.2 protein - Age (adjusted for sex and generation)
pval_protein_age <- ggplot(output, aes(x=p.Prot_Age.Sex)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="Protein/Age",
         x="p-value")
# 2.1 mRNA - Sex (adjusted for age and gernation)
pval_mRNA_sex <- ggplot(output, aes(x=p.mRNA_Sex.Age)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="mRNA/Sex",
         x="p-value")
# 2.2 protein - Sex (adjusted for age and gernation)
pval_protein_sex <- ggplot(output, aes(x=p.Prot_Sex.Age)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="Protein/Sex",
         x="p-value")
# 3.1 mRNA - Age:Sex interaction (age as variate)
pval_mRNA_int <- ggplot(output, aes(x=p.mRNA_Interaction)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="mRNA/(Age*Sex)",
         x="p-value")
# 3.2 protein - Age:Sex interaction (age as variate)
pval_protein_int <- ggplot(output, aes(x=p.Prot_Interaction)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="Protein/(Age*Sex)",
         x="p-value")

# Plot fdr
# 1.1 mRNA - Age (adjusted for sex and generation)
fdr_mRNA_age <- ggplot(output, aes(x=fdr.mRNA_Age.Sex)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="mRNA/Age",
         x="fdr")
# 1.2 protein - Age (adjusted for sex and generation)
fdr_protein_age <- ggplot(output, aes(x=fdr.Prot_Age.Sex)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="Protein/Age",
         x="fdr")
# 2.1 mRNA - Sex (adjusted for age and gernation)
fdr_mRNA_sex <- ggplot(output, aes(x=fdr.mRNA_Sex.Age)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="mRNA/Sex",
         x="fdr")
# 2.2 protein - Sex (adjusted for age and gernation)
fdr_protein_sex <- ggplot(output, aes(x=fdr.Prot_Sex.Age)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="Protein/Sex",
         x="fdr")
# 3.1 mRNA - Age:Sex interaction (age as variate)
fdr_mRNA_int <- ggplot(output, aes(x=fdr.mRNA_Interaction)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="mRNA/(Age*Sex)",
         x="fdr")
# 3.2 protein - Age:Sex interaction (age as variate)
fdr_protein_int <- ggplot(output, aes(x=fdr.Prot_Interaction)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="Protein/(Age*Sex)",
         x="fdr")

pdf("./Plot/mRNA-histograms.pdf", width = 12, height = 6)
grid.arrange(pval_mRNA_age, pval_mRNA_sex, pval_mRNA_int,
             fdr_mRNA_age, fdr_mRNA_sex, fdr_mRNA_int,
             ncol = 3)
dev.off()

pdf("./Plot/protein-histograms.pdf", width = 12, height = 6)
grid.arrange(pval_protein_age, pval_protein_sex, pval_protein_int,
             fdr_protein_age, fdr_protein_sex, fdr_protein_int,
             ncol = 3)
dev.off()
