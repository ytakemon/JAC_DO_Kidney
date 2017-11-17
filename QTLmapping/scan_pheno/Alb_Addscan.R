library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO1045_kidney.Rdata")

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, MM_snps, pos_column = "marker")
snps$chr <- as.character(MM_snps$chr)
MM_snps$chr[MM_snps$chr=="X"] <- "20"
map <- map_df_to_list(map = MM_snps, pos_column = "bp")
