plist <- as.numeric(commandArgs(trailingOnly = TRUE))

print(Sys.time())
cat(paste("Mapping", length(plist), "genes", "\n"))
print(plist)

library(qtl2)
library(qtl2convert)
library(tidyverse)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")

# subset to IntAgeChr == 12 list
Chr12_list <- readr::read_csv("./QTLscan/scanBestMarker_mrna/BestMarker_BestperGene_mrna_thr8.csv", guess_max = 4200) %>%
  filter(IntAgeChr == "12")
sub_annot.mrna <- annot.mrna %>% filter(id %in% Chr12_list$id) #dim(sub_annot.mrna) 153
sub_expr.mrna <- expr.mrna[,Chr12_list$id]

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")

# add Akt1 mRNA to annot.samples
annot.samples$Med <- expr.mrna[,annot.mrna[annot.mrna$symbol == "Akt1",]$id]

addcovar <- model.matrix(~ Sex + Age + Generation + Med, data=annot.samples)
intcovar <- model.matrix(~ Age, data = annot.samples)

plist <- plist[plist<=ncol(sub_expr.mrna)]
output <- sub_annot.mrna[plist,]
Max <- NULL
for (p in plist) {

  cat("Scanning ",which(p==plist)," out of ",length(plist),"\n")

  # scan additive model:
  LODadd <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=sub_expr.mrna[,p],
               addcovar=addcovar[,-1],
               cores=10,
               reml=TRUE)

  # scan full model:
  LODfull <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=sub_expr.mrna[,p],
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
    diff$IntAgeChr <- str_split_fixed(rownames(diff),"_",2)[,1] # get chr
    diff <- diff[diff$IntAgeChr %in% "12",]
    max <- diff[which(diff$IntAgeLODDiff == max(diff$IntAgeLODDiff, na.rm = TRUE)[1])[1],]
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

filename <- paste0("./QTLscan/scanBestMarker_mrna/maxLODscan_Akt1_m/maxLODscan_batch_",plist[1],".csv")
write_csv(output, path = filename)
print(Sys.time())
