# Identify animals that have the average given expression values for the given group.
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggsci)
library(stringr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")

# RNA-seq
identical(colnames(raw.mrna),annot.mrna$id)
mrna <- annot.mrna
mrna$expr <- NA
for (i in 1:nrow(annot.mrna)){
  mrna$expr[i] <- !all(raw.mrna[,i] < 1)
}
expressed <- mrna[mrna$duplicated == FALSE,]
expressed <- expressed[expressed$expr == TRUE,]
mrna <- raw.mrna
mrna <- mrna[,expressed$id]
write.csv(mrna, file = "./SuppData/RawExpressedRNA.csv", row.names = TRUE, quote = FALSE)
write.csv(expressed, file = "./SuppData/AnnotExpressedRNA.csv", row.names = FALSE, quote = FALSE)

# Proteomics
identical(colnames(raw.protein),annot.protein$id)
prot <- annot.protein
prot$expr <- NA
for (i in 1:nrow(annot.protein)){
  prot$expr[i] <- !all(raw.protein[,i] < 1)
}
expressed <- prot[prot$expr == TRUE,]
prot <- raw.protein
prot <- prot[, expressed$id]
write.csv(prot, file = "./SuppData/RawExpressedProt.csv", row.names = TRUE, quote = FALSE)
write.csv(expressed, file = "./SuppData/AnnotExpressedProt.csv", row.names = FALSE, quote = FALSE)

# RNA-seq : Proteomics Common
mrna <- annot.mrna
mrna$expr <- NA
for (i in 1:nrow(annot.mrna)){
  mrna$expr[i] <- !all(raw.mrna[,i] < 1)
}
mrna <- mrna[mrna$duplicated == FALSE,]
mrna <- mrna[mrna$expr == TRUE,]

prot <- annot.protein
prot$expr <- NA
for (i in 1:nrow(annot.protein)){
  prot$expr[i] <- !all(raw.protein[,i] < 1)
}
prot <- prot[prot$expr == TRUE,]

common <- prot[prot$gene_id %in% mrna$id,]
write.csv(common, file = "./SuppData/AnnotExpressedCommon.csv", row.names = FALSE, quote = FALSE)
