plist <- as.numeric(commandArgs(trailingOnly = TRUE))

print(Sys.time())
cat(paste("Mapping", length(plist), "genes", "\n"))
print(plist)

library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(tidyverse)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")

addcovar <- model.matrix(~ Sex + Age + Generation + Protein.Batch + Protein.Channel, data=annot.samples)
intcovar <- model.matrix(~ Age, data = annot.samples)

plist <- plist[plist<=ncol(expr.protein)]
output <- annot.protein[plist,]

MaxChr <- NULL
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
  maxChr <- function(lod){
    lod <- as.data.frame(lod)
    lod$chr <- str_split_fixed(rownames(lod),"_",2)[,1]
    lod$marker <- rownames(lod)
    maxChr <- lod %>% group_by(chr) %>% summarise(
        maxLOD = max(pheno1),
        maxMarker = marker[which(pheno1 == max(pheno1)[1])[1]]) %>% arrange(chr)
    return(maxChr)
  }

  MaxChr_add <- maxChr(LODadd)
  MaxChr_full <- maxChr(LODfull)

  DiffMaxChr <- MaxChr_full %>%
    rename(
      fullmaxLOD = maxLOD,
      fullmax_marker = maxMarker) %>%
    mutate(
      addmaxLOD = MaxChr_add$maxLOD,
      maxLODDiff = fullmaxLOD-addmaxLOD,
      fullmax_pos = NA)

  # Get max diff chr and annotate marker and pos
  maxchr <- DiffMaxChr %>% filter(maxLODDiff == max(maxLODDiff)[1])
  maxchr$fullmax_pos <- snps[snps$marker == maxchr$fullmax_marker,]$bp
  MaxChr <- rbind(MaxChr,maxchr) # } #test loop

  # save lod object
  filename_add <- paste0("./QTLscan/scanBestChr_protein/addscan/", annot.protein$id[p], "_", annot.protein$symbol[p], ".rds")
  filename_full <- paste0("./QTLscan/scanBestChr_protein/fullscan/", annot.protein$id[p], "_", annot.protein$symbol[p], ".rds")
  saveRDS(LODadd, file=filename_add)
  saveRDS(LODfull, file=filename_full)
}

output <- output %>% mutate(
  FullMaxLOD = MaxChr$fullmaxLOD,
  AddMaxLOD = MaxChr$addmaxLOD,
  IntAgeMaxChr = MaxChr$chr,
  IntAgeMaxPos = MaxChr$fullmax_pos,
  IntAgeMaxLODDiff = MaxChr$maxLODDiff
)

filename_intmax <- paste0("./QTLscan/scanBestChr_protein/intmaxchr/maxLODscan_batch_",plist[1],".csv")
write_csv(output, path = filename_intmax)
print(Sys.time())
