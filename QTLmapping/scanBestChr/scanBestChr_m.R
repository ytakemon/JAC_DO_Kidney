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

addcovar <- model.matrix(~ Sex + Age + Generation, data=annot.samples)
intcovar <- model.matrix(~ Age, data = annot.samples)

plist <- plist[plist<=ncol(expr.mrna)]
output <- annot.mrna[plist,]

MaxChr <- NULL
for (p in plist) {

  cat("Scanning ",which(p==plist)," out of ",length(plist),"\n")

  # scan additive model:
  LODadd <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=expr.mrna[,p],
               addcovar=addcovar[,-1],
               cores=10,
               reml=TRUE)

  # scan full model:
  LODfull <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=expr.mrna[,p],
               addcovar=addcovar[,-1],
               intcovar=intcovar[,-1],
               cores=10,
               reml=TRUE)

  # get max per chromosome:
  maxChr <- function(lod){
    lod <- as.data.frame(lod)
    lod$chr <- str_split_fixed(rownames(lod),"_",2)[,1]
    maxChr <- lod %>% group_by(chr) %>% summarise(max = max(pheno1)) %>% arrange(chr)
    return(maxChr)
  }

  MaxChr_add <- maxChr(LODadd)
  MaxChr_int <- maxChr(LODfull)

  DiffMaxChr <- MaxChr_int %>% rename(fullmax = max) %>%
    mutate(
      addmax = MaxChr_add$max,
      diffmax = fullmax-addmax
    )

  maxchr <- DiffMaxChr %>% filter(diffmax == max(diffmax))
  MaxChr <- rbind(MaxChr,maxchr)

  # save lod object
  filename_add <- paste0("./QTLscan/scanBestChr_mrna/addscan/", annot.mrna$id[p], "_", annot.mrna$symbol[p], ".rds")
  filename_full <- paste0("./QTLscan/scanBestChr_mrna/fullscan/", annot.mrna$id[p], "_", annot.mrna$symbol[p], ".rds")
  saveRDS(LODadd, file=filename_add)
  saveRDS(LODfull, file=filename_full)
}

output <- output %>% mutate(
  IntAgeMaxChr = MaxChr$chr,
  FullLODMax = MaxChr$fullmax,
  AddLODMax = MaxChr$addmax,
  IntAgeLODMaxDiff = MaxChr$diffmax
)

filename_intmax <- paste0("./QTLscan/scanBestChr_mrna/intmaxchr/maxLODscan_batch_",plist[1],".csv")
write_csv(output, path = filename_intmax)
print(Sys.time())
