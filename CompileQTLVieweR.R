# R/3.4.4
# Refer to: https://github.com/churchill-lab/qtl-viewer/blob/master/docs/QTLViewerDataStructures.md
library(qtl2)
library(qtl2convert)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney.RData")

# Need the following in this .Rdata file:
# ensembl.version
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

ensembl.version <- 82
genoprobs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
K <- calc_kinship(genoprobs, type = "loco")
map <- map_df_to_list(map = snps, pos_column = "bp")
markers <- snps
markers$pos <- markers$pos * 1e-6 # convert to megabases

# dataset.mrna: --------------------------------------------------------------
datatype <- "mRNA"

annots.mrna <- annot.mrna %>%
  filter(!(chr %in% c("MT","Y"))) %>%
  mutate(nearest.marker.id = NA)

for (i in 1:length(annots.mrna$nearest_snp)){
  annots.mrna$nearest.marker.id[i] <- markers$marker[annots.mrna$nearest_snp[i]]
}
annots.mrna <- annots.mrna %>%
  mutate(
    start = start * 1e-6, #convert to megabases
    end = end * 1e-6, #convert to megabases
    middle = middle_point * 1e-6, #convert to megabases
  ) %>%
  filter(duplicated == FALSE) %>%
  rename( gene_id = id) %>%
  select(-c(nearest_snp, middle_point))

rownames(annots.mrna) <- annots.mrna$gene_id


#expr <- expr.mrna
expr.mrna <- expr.mrna[,annots.mrna$gene_id]

# Find nearest marker for additive QTL lod peaks
additive <- read.csv("./QTLscan/output/eQTLAllAdditive.csv", header = TRUE, stringsAsFactors = FALSE) # additive only for now (from Matt)
colnames(additive)[1] <- "annot.id"
colnames(additive)[9] <- "lod"
additive$marker.id <- NA
additive$AdditivePos <- additive$AdditivePos * 1e-6 # Convert to megabases
# match nearest marker
for (i in 1:length(additive$marker.id)){
  chr <- additive$AdditiveChr[i]
  sub <- markers[markers$chr == chr,]
  additive$marker.id[i] <- sub$marker[which.min(abs(sub$pos - additive$AdditivePos[i]))]
}
#Check
is.data.frame(additive)

# interactive lod peaks
age_int <- read.csv("./QTLscan/output/eQTLAllInteractiveAge.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(age_int)[1] <- "annot.id"
colnames(age_int)[9] <- "lod"
age_int$marker.id <- NA
age_int$IntAgePos <- age_int$IntAgePos * 1e-6 # Convert to megabases
# match nearest marker
for (i in 1:length(age_int$marker.id)){
  chr <- age_int$IntAgeChr[i]
  sub <- markers[markers$chr == chr,]
  age_int$marker.id[i] <- sub$marker[which.min(abs(sub$pos - age_int$IntAgePos[i]))]
}
# Check
is.data.frame(age_int)
lod.peaks <- list("additive" = additive,
                  "age_int" = age_int)

samples <- as.data.frame(annot.samples)
colnames(samples)[1] <- "id"

# Grab raw data from EMASE (total gene expression): sample id x mrna annot.
#raw <- data.frame(matrix(NA, nrow = length(samples$id), ncol = nrow(annots.mrna)))
#rownames(raw) <- samples$id
#colnames(raw) <- rownames(annots.mrna)[order(rownames(annots.mrna))]
#emase_dir <- "/hpcdata/gac/projects/JAC_DO_Kidney_RNASeq/emase_m4_gbrs/"
#
#for (i in 1:nrow(raw)){
#  file <- paste0(emase_dir, rownames(raw)[i], ".emase.genes.effective_read_counts")
#  file <- read.delim(file, stringsAsFactors = FALSE)
#  file <- file[file$locus %in% rownames(annots.mrna),]
#
#  raw[i,] <- file$total
#  print(paste0(i," out of ", nrow(raw)))
#}
#raw <- as.matrix(raw)

# covars
covar <- covar[,-1]
covar.factors <- data.frame(
  column.name = c("Sex", "Age", "Generation"),
  display.name = c("Sex", "Age", "Generation"),
  int.covar = c(NA, "factor", NA),
  lod.peaks = c(NA, "age_int", NA),
  covar.name = c(NA,"Age", NA),
  stringsAsFactors = FALSE)


dataset.mrna <- list("annots" = annots.mrna,
                    "covar" = covar,
                    "covar.factors" = covar.factors,
                    "datatype" = datatype,
                    "rankz" = expr.mrna,
                    "lod.peaks" = lod.peaks,
                    "samples" = samples)

# save(genome.build, genoprobs, K, map, markers, dataset.mrna, file = "./RNAseq_data/DO188b_kidney_201711.Rdata")
# dataset.protein: -------------------------------------------------------------
datatype <- "protein"

annots.protein <- annot.protein %>%
  mutate(nearest.marker.id = NA)

for (i in 1:length(annots.protein$nearest_snp)){
  annots.protein$nearest.marker.id[i] <- markers$marker[annots.protein$nearest_snp[i]]
}
annots.protein <- annots.protein %>% mutate(
  start = start * 1e-6, # convert to megabases
  end = end * 1e-6, # convert to megabases
  middle = middle_point * 1e-6) %>%
  select(-c(nearest_snp, middle_point)) %>%
  rename(protein_id = id)
rownames(annots.protein) <- annots.protein$protein_id

# expression data
expr <- expr.protein

# Find nearest marker for additive QTL lod peaks
additive <- read.csv("./QTLscan/output/pQTLAllAdditive.csv", header = TRUE, stringsAsFactors = FALSE) # additive only for now (from Matt)
colnames(additive)[1] <- "annot.id"
colnames(additive)[9] <- "lod"
additive$marker.id <- NA
additive$AdditivePos <- additive$AdditivePos * 1e-6 # Convert to megabases
# match nearest marker
for (i in 1:length(additive$marker.id)){
  chr <- additive$AdditiveChr[i]
  sub <- markers[markers$chr == chr,]
  additive$marker.id[i] <- sub$marker[which.min(abs(sub$pos - additive$AdditivePos[i]))]
}
#Check
is.data.frame(additive)
# interactive lod peaks
age_int <- read.csv("./QTLscan/output/pQTLAllInteractiveAge.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(age_int)[1] <- "annot.id"
colnames(age_int)[9] <- "lod"
age_int$marker.id <- NA
age_int$IntAgePos <- age_int$IntAgePos * 1e-6 # Convert to megabases
# match nearest marker
for (i in 1:length(age_int$marker.id)){
  chr <- age_int$IntAgeChr[i]
  sub <- markers[markers$chr == chr,]
  age_int$marker.id[i] <- sub$marker[which.min(abs(sub$pos - age_int$IntAgePos[i]))]
}
# Check
is.data.frame(age_int)
lod.peaks <- list("additive" = additive,
                  "age_int" = age_int)

# covar
covar <- covar[,-1]
covar.factors <- data.frame(
  column.name = c("Sex", "Age", "Generation"),
  display.name = c("Sex", "Age", "Generation"),
  int.covar = c(NA, "factor", NA),
  lod.peaks = c(NA, "age_int", NA),
  covar.name = c(NA, "Age", NA),
  stringsAsFactors = FALSE)

dataset.protein <- list("annots" = annots.protein,
                     "covar" = covar,
                     "covar.factors" = covar.factors,
                     "datatype" = datatype,
                     "rankz" = expr.protein,
                     "lod.peaks" = lod.peaks,
                     "samples" = samples)

save(ensembl.version, genoprobs, K, map, markers, dataset.mrna, dataset.protein, file = "./RNAseq_data/DO188b_kidney_20180926_YT.Rdata")

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
