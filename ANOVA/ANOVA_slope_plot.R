# R/3.4.1
library(ggplot2)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("RNAseq_data/DO188b_kidney_noprobs.RData")
df <- read.csv("Anova_output/kidney_anova_slope_output.csv", header = T)
sig <- read.csv("Anova_output/kidney_anova_sig_output.csv", header = T)

mcols <- grep("^m.", colnames(df))
df <- df[, c(1:8, mcols)]

sigcols <- grep("p.", colnames(sig))
sig <- sig[, c(1:8,sigcols)]
# Plot age realted mRNA/Protein slopes
ggplot(df, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex)) +
      geom_point()
