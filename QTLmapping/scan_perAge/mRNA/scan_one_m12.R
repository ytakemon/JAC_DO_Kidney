plist <- as.numeric(commandArgs(trailingOnly = TRUE))

cat(paste("Mapping", length(plist), "genes", "\n"))
print(plist)

library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")

#subset to age group
genoprobs <- genoprobs[annot.samples$Age == "12",,]
expr.mrna <- expr.mrna[annot.samples$Age == "12",]
annot.samples <- annot.samples[annot.samples$Age == "12",]


# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")
K <- calc_kinship(probs, type = "loco", cores = 10)

plist <- plist[plist<=ncol(expr.mrna)]

for (p in plist) {

  cat("Scanning ",which(p==plist)," out of ",length(plist),"\n")

  addcovar <- model.matrix(~ Sex + Generation, data=annot.samples)

  # lod score
  lod <- scan1(genoprobs=probs,
               kinship=K,
               pheno=expr.mrna[,p],
               addcovar=addcovar[,-1], cores=10, reml=TRUE)

  # save lod object
  file_name <- paste0("./QTLscan/addscan_mrna_perAge/12mo/", annot.mrna$id[p], "_", annot.mrna$symbol[p], ".rds")
  saveRDS(lod, file=file_name)
}
