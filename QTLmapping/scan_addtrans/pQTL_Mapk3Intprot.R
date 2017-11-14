# run with:
# qsub -v script=pQTL_Mapk3Intprot Rsubmit_args.sh

library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO188b_kidney.RData")

# Get list of genes with trans eQTL
list <- read.csv("./QTLscan/output/Threshold8_pQTL_intAge_pbatch.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == 7, ]
list <- list$id # ENS protein id

# Identify gene name,
# Mapk3 gene expression is thought to mediate genes that have pQTL on Chr7
trans <- "Mapk3"
trans <- annot.protein[annot.protein$symbol == trans,]

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")

for (p in 1:length(list)){

  cat("Scanning ",p," out of ",length(list),"\n")
  addcovar <- model.matrix(~ Sex + Age + Generation + Protein.Batch + Protein.Channel + expr.protein[,trans$id], data=annot.samples)
  intcovar <- model.matrix(~ Age, data=annot.samples)

  p <- list[p]
  p <- annot.protein[annot.protein$id == p,]

  # Lod score
  lod <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=expr.protein[, p$id],
               addcovar=addcovar[,-1],
               intcovar=intcovar[,-1],
               cores=10, reml=TRUE)

  # Save lod object
  file_name <- paste0("./QTLscan/intscan_prot_Mapk3prot/", p$id, "_", p$symbol, ".rds")
  saveRDS(lod, file = file_name)
}
