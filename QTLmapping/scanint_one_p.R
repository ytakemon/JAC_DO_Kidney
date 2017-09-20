plist <- as.numeric(commandArgs(trailingOnly = TRUE))

cat(paste("Mapping", length(plist), "genes", "\n"))
print(plist)

library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
load("DO188b_kidney.RData")

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
    intcovar <- model.matrix(~ Sex, data=annot.samples[present, ])
    intcovar2<- model.matrix(~ Age, data=annot.samples[present, ])
  } else {
    addcovar <- model.matrix(~ Sex + Age, data=annot.samples[present, ])
    intcovar <- model.matrix(~ Sex, data=annot.samples[present, ])
    intcovar2<- model.matrix(~ Age, data=annot.samples[present, ])
  }  

  # lod score
  lod1 <- scan1(genoprobs=probs[present,],
               kinship=Glist2,
               pheno=expr.protein[present,p],
               addcovar=addcovar[,-1],
               intcovar=intcovar[,-1], 
               cores=10, reml=TRUE)
  # save lod object
  lod2 <- scan1(genoprobs=probs[present,],
               kinship=Glist2,
               pheno=expr.protein[present,p],
               addcovar=addcovar[,-1],
               intcovar=intcovar2[,-1], 
               cores=10, reml=TRUE)
  # save lod object
  file_name1 <- paste0("intscan_prot/Sex/", annot.protein$id[p], "_", annot.protein$symbol[p], ".rds")
  file_name2 <- paste0("intscan_prot/Age/", annot.protein$id[p], "_", annot.protein$symbol[p], ".rds")
  saveRDS(lod1, file=file_name1)
  saveRDS(lod2, file=file_name2)
}

