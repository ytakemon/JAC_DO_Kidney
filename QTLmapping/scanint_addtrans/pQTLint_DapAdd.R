# run with:
# qsub -v script=pQTLint_DapAdd Rsubmit_args.sh

library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO188b_kidney.RData")

# Get list of genes with trans eQTL
list <- read.csv("./QTLscan/output/Threshold8_pQTL_intAge.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == 15, ]
list <- list$symbol

# Identify gene name
# Dap was not detected in protein data, but RNA expression data is available
# will use RNA expression for now since there is a high possibilitiy that it could
# be a part of Quad I and III. (increase protein == increase mrna)
trans <- "Dap"
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
trans <- other.ids(trans, "mRNA")

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")

for (p in 1:length(list)){

  cat("Scanning ",p," out of ",length(list),"\n")
  addcovar <- model.matrix(~ Sex + Age + Generation + expr.mrna[,trans$id], data=annot.samples)
  intcovar <- model.matrix(~ Sex, data=annot.samples)

  p <- list[p]
  p <- other.ids(p, "protein")

  # Lod score
  lod <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=expr.protein[, p$protein_id],
               addcovar=addcovar[,-1],
               intcovar=intcovar[,-1],
               cores=10, reml=TRUE)

  # Save lod object
  file_name <- paste0("./QTLscan/intscan_prot_addDap/", p$id, "_", p$symbol, ".rds")
  saveRDS(lod, file = file_name)
}
