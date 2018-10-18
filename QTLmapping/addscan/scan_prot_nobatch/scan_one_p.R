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


plist <- plist[plist<=ncol(expr.protein)]

for (p in plist) {

  cat("Scanning ",which(p==plist)," out of ",length(plist),"\n")

  present <- !is.na(expr.protein[,p])
  for (j in 1:20)
    Glist2 <- Glist[[j]][present,present]
  if (length(unique(annot.samples$Generation[present]))>1){
    addcovar <- model.matrix(~ Sex + Age + Generation, data=annot.samples[present, ])
  } else {
    addcovar <- model.matrix(~ Sex + Age + Protein.Batch + Protein.Channel, data=annot.samples[present, ])
  }

  # lod score
  lod <- scan1(genoprobs=probs[present,],
               kinship=Glist2,
               pheno=expr.protein[present,p],
               addcovar=addcovar[,-1],
               cores=10, reml=TRUE)

  # save lod object
  file_name <- paste0("./QTLscan/addscan_prot_nobatch/", annot.protein$id[p], "_", annot.protein$symbol[p], ".rds")
  saveRDS(lod, file=file_name)
}
