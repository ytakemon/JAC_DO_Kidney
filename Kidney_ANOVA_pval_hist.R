# ANOVA for both RNA and protein:
# 1. Age (adjusted for sex and generation)
# 2. Sex (adjusted for age and generation)
# 3. Age:Sex interaction with age as a variate (adjust for generation).
# (Compute p-values and FDR using BH method)

# References:
# https://github.com/simecek/TheAgingProteome/blob/master/analysis/pval-histograms.Rmd

# Usage on Cadillac
# cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts
# qsub -v I=kidney_anova_output.csv,script=Kidney_ANOVA_pval_hist Rsubmit_args.R

# R/3.3.2
library(ggplot2)
library(gridExtra)

# Load data
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
args <- commandArgs(trailingOnly = TRUE)
output <- list.files(pattern = args[[1]], recursive = TRUE)
output <- read.csv(output[[1]])

# 1.1 mRNA - Age (adjusted for sex and generation)
mRNA_age <- ggplot(output, aes(x=p.mRNA_Age.Sex)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="mRNA/Age",
         x="p-value")
# 1.2 protein - Age (adjusted for sex and generation)
protein_age <- ggplot(output, aes(x=p.Prot_Age.Sex)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="Protein/Age",
         x="p-value")
# 2.1 mRNA - Sex (adjusted for age and gernation)
mRNA_sex <- ggplot(output, aes(x=p.mRNA_Sex.Age)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="mRNA/Sex",
         x="p-value")
# 2.2 protein - Sex (adjusted for age and gernation)
protein_sex <- ggplot(output, aes(x=p.Prot_Sex.Age)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="Protein/Sex",
         x="p-value")
# 3.1 mRNA - Age:Sex interaction (age as variate)
mRNA_int <- ggplot(output, aes(x=p.mRNA_Interaction)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="mRNA/(Age*Sex)",
         x="p-value")
# 3.2 protein - Age:Sex interaction (age as variate)
protein_int <- ggplot(output, aes(x=p.Prot_Interaction)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="Protein/(Age*Sex)",
         x="p-value")


pdf("./Plot/pval-histograms.pdf", width = 12, height = 6)
grid.arrange(mRNA_age, mRNA_sex, mRNA_int, protein_age, protein_sex, protein_int, ncol = 3)
dev.off()
