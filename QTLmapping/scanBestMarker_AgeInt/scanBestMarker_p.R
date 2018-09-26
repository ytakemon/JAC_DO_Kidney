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

plist <- plist[plist<=ncol(expr.protein)]
output <- annot.protein[plist,]

Max <- NULL
for (p in plist) {

  cat("Scanning ",which(p==plist)," out of ",length(plist),"\n")

  # scan additive model:
  LODadd <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=expr.protein[,p],
               addcovar=addcovar[,-1],
               cores=10,
               reml=TRUE)

  # scan full model:
  LODfull <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=expr.protein[,p],
               addcovar=addcovar[,-1],
               intcovar=intcovar[,-1],
               cores=10,
               reml=TRUE)

  # get max per chromosome:
  maxMarker <- function(full, add){
    diff <- as.data.frame(full)
    colnames(diff) <- "FullLOD"
    diff$AddLOD <- add[,1]
    diff$IntAgeLODDiff <- diff$FullLOD - diff$AddLOD
    max <- diff[which(diff$IntAgeLODDiff == max(diff$IntAgeLODDiff, na.rm = TRUE)[1])[1],]
    max$IntAgeChr <- str_split_fixed(rownames(max),"_",2)[,1]
    max$IntAgePos <- snps[snps$marker == rownames(max),]$bp
    return(max)
  }
  max <- maxMarker(LODfull, LODadd)
  Max <- rbind(Max,max) # } #test loop
}

output <- output %>% mutate(
  FullLOD = Max$FullLOD,
  AddLOD = Max$AddLOD,
  IntAgeChr = Max$IntAgeChr,
  IntAgePos = Max$IntAgePos,
  IntAgeLODDiff = Max$IntAgeLODDiff
)

filename <- paste0("./QTLscan/scanBestMarker_protein/maxLODscan/maxLODscan_batch_",plist[1],".csv")
write_csv(output, path = filename)
print(Sys.time())
