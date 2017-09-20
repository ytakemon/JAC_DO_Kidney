plist <- as.numeric(commandArgs(trailingOnly = TRUE))

cat(paste("Mapping", length(plist), "genes", "\n"))
print(plist)

library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
load("DO188_kidney.RData")

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
map <- map_df_to_list(map = snps, pos_column = "bp")



plist <- plist[plist<=ncol(expr.protein)]

for (p in plist) {

  cat("Scanning ",which(p==plist)," out of ",length(plist),"\n")  

  # lod score
  lod <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=expr.protein[,p],
               addcovar=covar[,-1], cores=1, reml=TRUE)
  
  # save lod object
  file_name <- paste0("scanone_protein/", annot.protein$id[p], "_", annot.protein$symbol[p], ".rds")
  saveRDS(lod, file=file_name)
}

