plist <- as.numeric(commandArgs(trailingOnly = TRUE))

print(Sys.time())
cat(paste("Mapping", length(plist), "genes", "\n"))
print(plist)

library(qtl2)
library(qtl2convert)
library(tidyverse)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")
WB <- read.delim("./Phenotype/phenotypes/JAC_WB_kidney_ERK.txt")
pheno <- "Phospho_ERK1"

# Cleanup data and subset to match samples
# match pheno to samples
WB_clean <- WB[complete.cases(WB[,pheno]),] %>%
  rename(Mouse.ID = ID) %>%
  mutate(duplicated = (duplicated(Mouse.ID) | duplicated(Mouse.ID, fromLast = TRUE))) %>%
  filter((Mouse.ID %in% annot.samples$Mouse.ID) & duplicated == FALSE) %>%
  arrange(Mouse.ID)

# subset to IntAgeChr == 7 list
Chr7_list <- readr::read_csv("./QTLscan/scanBestMarker_protein/BestMarker_BestperGene_protein_thr8.csv", guess_max = 4200) %>%
  filter(IntAgeChr == "7")

sub_annot.protein <- filter(annot.protein, id %in% Chr7_list$id) #dim(sub_annot.protein) 262 genes
sub_expr.protein <- expr.protein[,Chr7_list$id]
sub_expr.protein <- sub_expr.protein[rownames(sub_expr.protein) %in% WB_clean$Mouse.ID,]
sub_annot.samples <- annot.samples[rownames(annot.samples) %in% WB_clean$Mouse.ID,]

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")

# add Akt1 mRNA to annot.samples
# add ERK1/Mapk3 mRNA to annot.samples
sub_annot.samples$Med <- WB_clean[,pheno]

addcovar <- model.matrix(~ Sex + Age + Generation + Med, data=sub_annot.samples)
intcovar <- model.matrix(~ Age, data = annot.samples)

plist <- plist[plist<=ncol(sub_expr.protein)]
output <- sub_annot.protein[plist,]
Max <- NULL
for (p in plist) {

  cat("Scanning ",which(p==plist)," out of ",length(plist),"\n")

  # scan additive model:
  LODadd <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=sub_expr.protein[,p],
               addcovar=addcovar[,-1],
               cores=10,
               reml=TRUE)

  # scan full model:
  LODfull <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=sub_expr.protein[,p],
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

filename <- paste0("./QTLscan/scanBestMarker_protein/maxLODscan_Erk1_WBtotal/maxLODscan_batch_",plist[1],".csv")
write_csv(output, path = filename)
print(Sys.time())
