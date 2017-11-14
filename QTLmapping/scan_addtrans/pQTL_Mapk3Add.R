# run with:
# qsub -v script=pQTL_add_batch Rsubmit_args.sh

library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO188b_kidney.RData")

# Get list of genes with trans eQTL
list <- read.csv("./QTLscan/output/pQTLAllInteractiveAge_pbatch.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == 7, ]
list <- list$symbol

# Identify gene name
# Mapk3 was not detected in protein data, but RNA expression data is available
# will use RNA expression for now since there is a high possibilitiy that it could
# be a part of Quad I and III. (increase protein == increase mrna)
trans <- "Mapk3"
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
trans <- other.ids(trans, "protein")


# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")

for (p in 1:length(list)){

  cat("Scanning ",p," out of ",length(list),"\n")
  addcovar <- model.matrix(~ Sex + Age + Generation + Protein.Batch + Protein.Channel + expr.mrna[,trans$id], data=annot.samples)

  p <- list[p]
  p <- other.ids(p, "protein")

  # Lod score
  lod <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=expr.protein[, p$protein_id],
               addcovar=addcovar[,-1],
               cores=10, reml=TRUE)

  # Save lod object
  file_name <- paste0("./QTLscan/addscan_prot_batch/", p$protein_id, "_", p$symbol, ".rds")
  saveRDS(lod, file = file_name)
}
