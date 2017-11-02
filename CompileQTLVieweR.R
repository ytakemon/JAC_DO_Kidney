# R/3.4.1
# Refer to: https://github.com/churchill-lab/qtl-viewer/blob/master/docs/QTLViewerDataStructures.md
# Compiling Rdata: ./DO188b_kidney_201711.Rdata
library(qtl2)
library(qtl2convert)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney.RData")
load("./RNAseq_data/DO188b_kidney_201711.Rdata")

# Need the following in this .Rdata file:
# genome.build
# genoprobs
# K
# map
# markers
# dataset.mrna:
#   annots
#   covar
#   covar.factors
#   datatype
#   display.name
#   ensembl.version
#   expr
#   lod.peaks
#   samples
# dataset.protein: (see dataset.mrna for list)
# dataset.phenotype: (see dataset.mrna for list)

genome.build <- as.character("GRCm38")
genoprobs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
K <- calc_kinship(genoprobs, type = "loco")
map <- map_df_to_list(map = snps, pos_column = "bp")
markers <- snps
markers$pos <- markers$pos * 1e-6 # convert to megabases

# dataset.mrna: --------------------------------------------------------------
datatype <- "mRNA"

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

covar <- covar[,-1]
covar.factors <- data.frame(column.name = colnames(covar),
                            display.name = c("Sex", "Age", "G11","G12","G8","G9","Sex*Age"))
expr <- expr.mrna

# Find nearest marker for additive QTL lod peaks
lod.peaks <- read.csv("./QTLscan/output/eQTLAllAdditive.csv", header = TRUE, stringsAsFactors = FALSE) # additive only for now (from Matt)
colnames(lod.peaks)[1] <- "annot.id"
lod.peaks$marker.id <- NA
lod.peaks$AdditivePos <- lod.peaks$AdditivePos * 1e-6 # Convert to megabases
# match nearest marker
for (i in 1:length(lod.peaks$marker.id)){
  chr <- lod.peaks$AdditiveChr[i]
  sub <- markers[markers$chr == chr,]
  lod.peaks$marker.id[i] <- sub$marker[which.min(abs(sub$pos - lod.peaks$AdditivePos[i]))]
}

raw <- matrix(NA, nrow = nrow(samples), ncol = nrow(annots.mrna),
              dimnames = list(c(rownames(samples)), c(rownames(annots.mrna))))




R data type: matrix

Description: raw expression data, samples (rows) x annots (columns)















samples <- annot.samples
colnames(samples)[1] <- "id"

dataset.yuka <- c(“annots”=annots.mrna, “covar”=covar…

# Saved so far (just in case)
# save(genome.build, genoprobs, K, map, markers, annots.mrna, file = "./RNAseq_data/DO188b_kidney_201711.Rdata")
