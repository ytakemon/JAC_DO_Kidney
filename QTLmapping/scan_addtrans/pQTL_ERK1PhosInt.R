# run with:
# qsub -v script=pQTL_ERK1PhosInt Rsubmit_args.sh

library(qtl2)
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO188b_kidney.RData")
erk1 <- read.delim("./Phenotype/phenotypes/JAC_WB_kidney_ERK.txt")
pheno <- "Phospho_ERK1"

# Cleanup data and subset to match samples
erk1 <- erk1[complete.cases(erk1[,pheno]),]
erk1 <- erk1[erk1$ID %in% annot.samples$Mouse.ID, ]
erk1$duplicated <- (duplicated(erk1$ID) | duplicated(erk1$ID, fromLast = TRUE))
erk1 <- erk1[erk1$duplicated == FALSE,]
erk1 <- erk1[erk1$ID %in% rownames(genoprobs),] # 55 animals (33 females, 22 males)
rownames(erk1) <- erk1$ID

# Get list of genes with trans eQTL
list <- read.csv("./QTLscan/output/Threshold6_pQTL_intAge_pbatch.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == 7, ]
list <- list$id # ENS protein id

# match animal data
genoprobs <- genoprobs[rownames(genoprobs) %in% rownames(erk1),,]
annot.samples <- annot.samples[rownames(annot.samples) %in% rownames(erk1),]
expr.mrna <- expr.mrna[rownames(expr.mrna) %in% rownames(erk1),]
stopifnot(identical(rownames(genoprobs), rownames(erk1)))

# add erk to annotation
annot.samples$erk1 <- erk1[,pheno]

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")

for (p in 1:length(list)){

  cat("Scanning ",p," out of ",length(list),"\n")
  addcovar <- model.matrix(~ Sex + Age + Generation + Protein.Batch + Protein.Channel + log(erk1), data=annot.samples)
  intcovar <- <- model.matrix(~ Age, data=annot.samples)
  p <- list[p]
  p <- annot.protein[annot.protein$id == p,]

  # Lod score
  lod <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=expr.protein[, p$id],
               addcovar=addcovar[,-1],
               intcovar=intcovar[,-1],
               cores=10, reml=TRUE)

  # Save lod object
  file_name <- paste0("./QTLscan/intscan_prot_ERK1Phos/", p$id, "_", p$symbol, ".rds")
  saveRDS(lod, file = file_name)
}
