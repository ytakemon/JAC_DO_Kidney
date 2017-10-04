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
args <- commandArgs(trailingOnly = TRUE) # args <- "kidney_anova_factor_table.csv"
output <- list.files(path = "./Anova_output/", pattern = paste0("^",args[[1]]), recursive = TRUE)
output <- read.csv(paste0("./Anova_output/",output[[1]]), header = T)

# Plot pvalues
# 1.1 mRNA - Age (adjusted for sex and generation) ----------------------------
pval_mRNA_Age6_12 <- ggplot(output, aes(x=p.mRNA_Age.Sex_fact_6_12)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="mRNA/Age 6mo - 12mo",
         x="p-value")
pval_mRNA_Age12_18 <- ggplot(output, aes(x=p.mRNA_Age.Sex_fact_12_18)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="mRNA/Age 12mo - 18mo",
         x="p-value")
pval_mRNA_Age6_18 <- ggplot(output, aes(x=p.mRNA_Age.Sex_fact_6_18)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="mRNA/Age 6mo - 18mo",
         x="p-value")

# 1.2 Protein - Age (adjusted for sex and generation)---------------------------
pval_Prot_Age6_12 <- ggplot(output, aes(x=p.Prot_Age.Sex_fact_6_12)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="Prot/Age 6mo - 12mo",
         x="p-value")
pval_Prot_Age12_18 <- ggplot(output, aes(x=p.Prot_Age.Sex_fact_12_18)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="Prot/Age 12mo - 18mo",
         x="p-value")
pval_Prot_Age6_18 <- ggplot(output, aes(x=p.Prot_Age.Sex_fact_6_18)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="Prot/Age 6mo - 18mo",
         x="p-value")
# 2.1 mRNA - Age*Sex (adjusted for sex and generation) -------------------------
pval_mRNA_Age_Int_6_12 <- ggplot(output, aes(x=p.mRNA_Int_fact_6_12)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="mRNA/Age*Sex 6mo - 12mo",
         x="p-value")
pval_mRNA_Age_Int_12_18 <- ggplot(output, aes(x=p.mRNA_Int_fact_12_18)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="mRNA/Age*Sex 12mo - 18mo",
         x="p-value")
pval_mRNA_Age_Int_6_18 <- ggplot(output, aes(x=p.mRNA_Int_fact_6_18)) +
    geom_histogram(binwidth=0.04) +
    theme_bw() +
    labs(title="Kidney",
         subtitle="mRNA/Age*Sex 6mo - 18mo",
         x="p-value")







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
