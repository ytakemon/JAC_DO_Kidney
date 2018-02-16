# R/3.4.1
# Refer to: https://github.com/churchill-lab/qtl-viewer/blob/master/docs/QTLViewerDataStructures.md
# Compiling Rdata: ./DO188b_kidney_201711.Rdata
library(qtl2)
library(qtl2convert)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney.RData")
# load("./RNAseq_data/DO188b_kidney_201711.Rdata")

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
#expr <- expr.mrna
expr.mrna <- expr.mrna[,annots.mrna$id]

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

samples <- as.data.frame(annot.samples)
colnames(samples)[1] <- "id"

# Grab raw data from EMASE (total gene expression): sample id x mrna annot.
raw <- data.frame(matrix(NA, nrow = length(samples$id), ncol = nrow(annots.mrna)))
rownames(raw) <- samples$id
colnames(raw) <- rownames(annots.mrna)[order(rownames(annots.mrna))]
emase_dir <- "/hpcdata/gac/projects/JAC_DO_Kidney_RNASeq/emase_m4_gbrs/"

for (i in 1:nrow(raw)){
  file <- paste0(emase_dir, rownames(raw)[i], ".emase.genes.effective_read_counts")
  file <- read.delim(file, stringsAsFactors = FALSE)
  file <- file[file$locus %in% rownames(annots.mrna),]

  raw[i,] <- file$total
  print(paste0(i," out of ", nrow(raw)))
}
raw <- as.matrix(raw)

dataset.mrna <- list("annots" = annots.mrna,
                    "covar" = covar,
                    "covar.factors" = covar.factors,
                    "datatype" = datatype,
                    "expr" = expr.mrna,
                    "lod.peaks" = lod.peaks,
                    "raw" = raw,
                    "samples" = samples)

# save(genome.build, genoprobs, K, map, markers, dataset.mrna, file = "./RNAseq_data/DO188b_kidney_201711.Rdata")
# dataset.protein: -------------------------------------------------------------
datatype <- "protein"
annots.protein <- annot.protein
annots.protein$nearest.marker.id <- NA
for (i in 1:length(annots.protein$nearest_snp)){
  annots.protein$nearest.marker.id[i] <- markers$marker[annots.protein$nearest_snp[i]]
}
annots.protein$start <- annots.protein$start * 1e-6 # convert to megabases
annots.protein$end <- annots.protein$end * 1e-6 # convert to megabases
annots.protein$middle <- annots.protein$middle_point * 1e-6 # convert to megabases
rownames(annots.protein) <- annots.protein$id
annots.protein <- annots.protein %>% select(-c(nearest_snp, middle_point))

expr <- expr.protein

# Find nearest marker for additive QTL lod peaks
lod.peaks <- read.csv("./QTLscan/output/pQTLAllAdditive.csv", header = TRUE, stringsAsFactors = FALSE) # additive only for now (from Matt)
colnames(lod.peaks)[1] <- "annot.id"
lod.peaks$marker.id <- NA
lod.peaks$AdditivePos <- lod.peaks$AdditivePos * 1e-6 # Convert to megabases
# match nearest marker
for (i in 1:length(lod.peaks$marker.id)){
  chr <- lod.peaks$AdditiveChr[i]
  sub <- markers[markers$chr == chr,]
  lod.peaks$marker.id[i] <- sub$marker[which.min(abs(sub$pos - lod.peaks$AdditivePos[i]))]
}

dataset.protein <- list("annots" = annots.protein,
                     "covar" = covar,
                     "covar.factors" = covar.factors,
                     "datatype" = datatype,
                     "expr" = expr.protein,
                     "lod.peaks" = lod.peaks,
                     "samples" = samples)

# save(genome.build, genoprobs, K, map, markers, dataset.mrna, dataset.protein, file = "./RNAseq_data/DO188b_kidney_201802_YT.Rdata")

# Check format -----------------------------------------------------------------
# Run check: https://github.com/churchill-lab/qtl-viewer/blob/master/scripts/qtlDataCheck.R
rm(list = ls())
load("./RNAseq_data/DO188b_kidney_201802_YT.Rdata")
source("./Scripts/qtlDataCheck.R")

CheckVariables()
CheckDatasets() # requirements are differnt
CheckExtraVars(ls())
CheckDataNames(ds = get(apropos("^dataset")[1]))

# > source("./Scripts/qtlDataCheck.R")
# > CheckVariables()
# > CheckDatasets()
# [1] "ERROR: covar is type: matrix, should be type: data.frame, in dataset: dataset.mrna"
# [1] "ERROR: display.name is type: NULL, should be type: character, in dataset: dataset.mrna"
# [1] "ERROR: ensembl.version is type: NULL, should be type: numeric, in dataset: dataset.mrna"
# [1] "ERROR: covar is type: matrix, should be type: data.frame, in dataset: dataset.protein"
# [1] "ERROR: display.name is type: NULL, should be type: character, in dataset: dataset.protein"
# [1] "ERROR: ensembl.version is type: NULL, should be type: numeric, in dataset: dataset.protein"
# [1] "ERROR: raw is type: NULL, should be type: matrix, in dataset: dataset.protein"
# > CheckExtraVars(ls())
# [1] "Warning: the following extra variables were found..."
# [1] "CheckDataNames" "CheckDatasets"  "CheckExtraVars" "CheckVariables"
# [5] "IsVariableOK"
# > CheckDataNames(ds = get(apropos("^dataset")[1]))
# NULL
