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
Chr7_list <- readr::read_csv("./QTLscan/scanBestMarker_protein/BestMarker_BestperGene_protein_thr8.csv", guess_max = 4200) %>%
  filter(IntAgeChr == "7") %>% arrange(id)
Chr7_list$markers <- paste0(Chr7_list$IntAgeChr, "_",Chr7_list$IntAgePos * 1e6)

# subset data
sub_annot.protein <- annot.protein %>% filter(id %in% Chr7_list$id) #dim(sub_annot.protein) 153
sub_expr.protein <- expr.protein[,Chr7_list$id]

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")

# add Akt1 mRNA to annot.samples
annot.samples$Med <- expr.mrna[,annot.protein[annot.protein$symbol == "Mapk3",]$gene_id]

addcovar <- model.matrix(~ Sex + Age + Generation + Med, data=annot.samples)
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
    diff$IntAgeChr <- str_split_fixed(rownames(diff),"_",2)[,1] # get chr
    max <- diff[rownames(diff) %in% Chr7_list$markers[p],] # Identify marker from original scan
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

filename <- paste0("./QTLscan/scanBestMarker_protein/maxLODscan_Erk1_m/maxLODscan_batch_",plist[1],".csv")
write_csv(output, path = filename)
print(Sys.time())
