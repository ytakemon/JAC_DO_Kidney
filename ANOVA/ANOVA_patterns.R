# One question in my mind is whether all of these RNA/proteins are changing in the same mice?
# Are there groups of genes that always change together?

# R.3.4.1
library(ggplot2)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")
sig_df <- read.csv("./Anova_output/kidney_anova_sig_output.csv", stringsAsFactors = FALSE)

# Identify mice that have signigicant changes in mRNA/ Protein with Age
# Identify sig genes
# ~ Age*Sex (interaction)
mRNA_int_anova <- sig_df[sig_df$sig.mRNA_Interaction == TRUE,] #n = 0
Prot_int_anova <- sig_df[sig_df$sig.Prot_Interaction == TRUE,] #n = 45
# ~ Age (adj for sex, generation)
mRNA_age_anova <- sig_df[sig_df$sig.mRNA_Age.Sex == TRUE,][, c("gene_id", "p.mRNA_Age.Sex", "r.mRNA_Age.Sex")] # n = 163
Prot_age_anova <- sig_df[sig_df$sig.Prot_Age.Sex == TRUE,][, c("gene_id", "p.Prot_Age.Sex", "r.Prot_Age.Sex")] # n = 3940
# Duplicates are created here because multiple isoforms will correspond to the same gene_id
# For now remove duplicate genes
mRNA_age_anova <- mRNA_age_anova[!duplicated(mRNA_age_anova),]
Prot_age_anova <- Prot_age_anova[!duplicated(Prot_age_anova),]

# For PANTHER GO processing ----------------------------------------------------
# separate sig increase and decerased genes using coef values
mRNA_age_anova_pos <- mRNA_age_anova[mRNA_age_anova$r.mRNA_Age.Sex > 0,]
mRNA_age_anova_neg <- mRNA_age_anova[mRNA_age_anova$r.mRNA_Age.Sex < 0,]

Prot_age_anova_pos <- Prot_age_anova[Prot_age_anova$r.Prot_Age.Sex > 0,]
Prot_age_anova_neg <- Prot_age_anova[Prot_age_anova$r.Prot_Age.Sex < 0,]

# Export for Panther GO annoation
write.table(mRNA_age_anova_pos$gene_id,
            file = "./Anova_output/gene_lists/mRNA_age_sig_pos.csv",
            sep =",",
            quote = FALSE,
            row.names = FALSE)
write.table(mRNA_age_anova_neg$gene_id,
            file = "./Anova_output/gene_lists/mRNA_age_sig_neg.csv",
            sep =",",
            quote = FALSE,
            row.names = FALSE)
write.table(Prot_age_anova_pos$gene_id,
            file = "./Anova_output/gene_lists/Prot_age_sig_pos.csv",
            sep =",",
            quote = FALSE,
            row.names = FALSE)
write.table(Prot_age_anova_neg$gene_id,
            file = "./Anova_output/gene_lists/Prot_age_sig_neg.csv",
            sep =",",
            quote = FALSE,
            row.names = FALSE)

# See if significant genes between mRNA and Protein overlap --------------------
# look up from longer list (protein)
common <- unique(Prot_age_anova$gene_id[Prot_age_anova$gene_id %in% mRNA_age_anova$gene_id])
# There are 107
# Are these going up or down with age in mRNA/Protein
common_mRNA_age_anova_pos <- mRNA_age_anova_pos[mRNA_age_anova_pos$gene_id %in% common,] #n = 74
common_mRNA_age_anova_neg <- mRNA_age_anova_neg[mRNA_age_anova_neg$gene_id %in% common,] #n = 33
common_Prot_age_anova_pos <- Prot_age_anova_pos[Prot_age_anova_pos$gene_id %in% common,] #n = 81 (duplicate!)
common_Prot_age_anova_neg <- Prot_age_anova_neg[Prot_age_anova_neg$gene_id %in% common,] #n = 27

# Find duplicate in protein:-------------------------------------------
## x <- Prot_age_anova_pos$gene_id[Prot_age_anova_pos$gene_id %in% common]
## y <- Prot_age_anova_neg$gene_id[Prot_age_anova_neg$gene_id %in% common]
## common_x <- common[common %in% x] # 80 <- duplicate?
## common_y <- common[common %in% y]
## duplicate <- common_Prot_age_anova_pos[common_Prot_age_anova_pos$gene_id == common_Prot_age_anova_pos$gene_id[duplicated(common_Prot_age_anova_pos$gene_id)],]
# duplicate found in protein data, has unique pval and coef value. Keep.

# Are sig postive mRNA gene always postive Protein? This should add up to 74
pos_m_p <- common_mRNA_age_anova_pos[common_mRNA_age_anova_pos$gene_id %in% common_Prot_age_anova_pos$gene_id,]
not_pos_m_p <- common_mRNA_age_anova_pos[!(common_mRNA_age_anova_pos$gene_id %in% common_Prot_age_anova_pos$gene_id),]
# Are sig neg mRNA gene always neg Protein? This should add up to 33
neg_m_p <- common_mRNA_age_anova_neg[common_mRNA_age_anova_neg$gene_id %in% common_Prot_age_anova_neg$gene_id,]
not_neg_m_p <- common_mRNA_age_anova_neg[!(common_mRNA_age_anova_neg$gene_id %in% common_Prot_age_anova_neg$gene_id),]

# Process genes through PANTHER
write.table(pos_m_p$gene_id,
            file = "./Anova_output/gene_lists/Increase_age.csv",
            sep =",",
            quote = FALSE,
            row.names = FALSE)
write.table(neg_m_p$gene_id,
            file = "./Anova_output/gene_lists/Decrease_age.csv",
            sep =",",
            quote = FALSE,
            row.names = FALSE)
