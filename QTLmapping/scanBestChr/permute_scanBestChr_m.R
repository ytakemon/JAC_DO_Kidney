plist <- as.numeric(commandArgs(trailingOnly = TRUE))

print(Sys.time())
cat(paste("Mapping", length(plist), "genes", "\n"))
print(plist)

library(qtl2)
library(qtl2convert)
library(tidyverse)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")

addcovar <- model.matrix(~ Sex + Age + Generation, data=annot.samples)
intcovar <- model.matrix(~ Age, data = annot.samples)

plist <- plist[plist<=ncol(expr.mrna)]
output_add <- output_full <- annot.mrna[plist,]
ADDPERM <- FULLPERM <- NULL
for (p in plist) {

  cat("Scanning ",which(p==plist)," out of ",length(plist),"\n")

  # permutation test
  addperm <- scan1perm(genoprobs = probs,
                      pheno=expr.mrna[,p],
                      kinship = Glist,
                      addcovar = addcovar[,-1],
                      n_perm = 10,
                      cores = 10,
                      reml = TRUE)

  fullperm <- scan1perm(genoprobs = probs,
                      pheno=expr.mrna[,p],
                      kinship = Glist,
                      addcovar = addcovar[,-1],
                      intcovar = intcovar[,-1],
                      n_perm = 10,
                      cores = 10,
                      reml = TRUE)

  # combine
  ADDPERM <- rbind(ADDPERM, addperm[,1])
  FULLPERM <- rbind(FULLPERM, fullperm[,1])
}

# add colnames and combine
colnames(ADDPERM) <- colnames(FULLPERM) <- c("perm1","perm2","perm3","perm4","perm5","perm6","perm7","perm8","perm9","perm10")
output_add <- cbind(output_add,ADDPERM)
output_full <- cbind(output_full,ADDPERM)

file_add <- paste0("./QTLscan/scanBestChr_mrna/perm/add/perm10scan_batch_",plist[1],".csv")
write_csv(output_add, path = file_add)

file_full <- paste0("./QTLscan/scanBestChr_mrna/perm/full/perm10scan_batch_",plist[1],".csv")
write_csv(output_full, path = file_full)
print(Sys.time())
