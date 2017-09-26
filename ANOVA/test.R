# One question in my mind is whether all of these RNA/proteins are changing in the same mice?
# Are there groups of genes that always change together?

# R.3.4.1
library(ggplot2)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")
sig_df <- read.csv("./Anova_output/kidney_anova_sig_output.csv", stringsAsFactors = FALSE)

# Identify mice that have signigicant changes in mRNA/ Protein with Age
# Identify sig genes
mRNA_age_anova <- sig_df[sig_df$sig.mRNA_Age.Sex == TRUE,][, c("gene_id", "p.mRNA_Age.Sex", "r.mRNA_Age.Sex")] # n = 166
Prot_age_anova <- sig_df[sig_df$sig.Prot_Age.Sex == TRUE,][, c("id", "p.Prot_Age.Sex", "r.Prot_Age.Sex")] # n = 3940
# pull out sig genes from animal expression data
expr_mRNA_sig <- expr_mRNA_sig[, mRNA_age_anova$gene_id]
expr_Prot_sig <- expr.protein[, Prot_age_anova$id]
# Add age info
expr_mRNA_sig <- cbind(annot.samples, expr_mRNA_sig)
expr_Prot_sig <- cbind(annot.samples, expr_Prot_sig)

#Example
# x <- 100
# ggplot(expr_mRNA_sig, aes(x = as.factor(Age), y = expr_mRNA_sig[,mRNA_age_anova$gene_id[x]], fill = as.factor(Age))) +
#        geom_boxplot(alpha = 0.5)+
#        geom_point()+
#        scale_x_discrete("Age")+
#        scale_y_continuous(paste0("Ranked ", mRNA_age_anova$gene_id[x], " mRNA Expression"))

# 
