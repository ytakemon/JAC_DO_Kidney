plist <- as.numeric(commandArgs(trailingOnly = TRUE))

cat(paste("Mapping", length(plist), "genes", "\n"))
print(plist)

library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")


plist <- plist[plist<=ncol(expr.mrna)]

for (p in plist) {

  cat("Scanning ",which(p==plist)," out of ",length(plist),"\n")

  addcovar <- model.matrix(~ Sex + Age + Generation, data=annot.samples)

  # lod score
  lod <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=expr.mrna[,p],
               addcovar=addcovar[,-1], cores=10, reml=TRUE)

  # save lod object
  file_name <- paste0("./QTLscan/addscan_mrna/", annot.mrna$id[p], "_", annot.mrna$symbol[p], ".rds")
  saveRDS(lod, file=file_name)
}
