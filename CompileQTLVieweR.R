# R/3.4.1
# Refer to: https://github.com/churchill-lab/qtl-viewer/blob/master/docs/QTLViewerDataStructures.md
# Compiling Rdata: ./DO188b_kidney_201711.Rdata
library(qtl2)
library(qtl2convert)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney.RData")

genome.build <- as.character("GRCm38")
genoprobs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
K <- calc_kinship(genoprobs, type = "loco")
map <- map_df_to_list(map = snps, pos_column = "bp")
markers <- snps
markers$pos <- markers$pos * 1e-6 # convert to megabases
# save(genome.build, genoprobs, K, map, markers, file = "./RNAseq_data/DO188b_kidney_201711.Rdata")

# dataset.mrna:
annots.mrna <- annot.mrna
annots.mrna <- annots.mrna[!(annots.mrna$chr %in% c("MT","Y")),]
annots.mrna$nearest.marker.id <- NA
for (i in 1:length(annots.mrna$nearest_snp)){
  annots.mrna$nearest.marker.id[i] <- markers$marker[annots.mrna$nearest_snp[i]]
}
annots.mrna$start <- annots.mrna$start * 1e-6 # convert to megabases
annots.mrna$end <- annots.mrna$end * 1e-6 # convert to megabases
annots.mrna$middle <- annots.mrna$middle_point * 1e-6 # convert to megabases
annots.mrna <- annots.mrna[annots.mrna$duplicated == FALSE,]
rownames(annots.mrna) <- annots.mrna$id
annots.mrna <- annots.mrna %>% select(-c(nearest_snp, middle_point))
# save(genome.build, genoprobs, K, map, markers, annots.mrna, file = "./RNAseq_data/DO188b_kidney_201711.Rdata")

covar <- covar[,-1]
