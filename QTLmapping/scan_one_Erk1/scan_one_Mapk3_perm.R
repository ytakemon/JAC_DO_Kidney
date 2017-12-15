# qsub -v script=scan_one_Mapk3_perm Rsubmit_args.sh

library(qtl2)
library(qtl2convert)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")
pheno <- "MAPK3_proteome"

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "pos")
K <- calc_kinship(probs, "loco", cores=3)
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="20"] <- "X"
map <- map_df_to_list(map = snps, pos_column = "pos")

addcovar <- model.matrix(~ Sex + Generation + Age + Protein.Batch + Protein.Channel , data = annot.samples)

prot <- "Mapk3"
prot <- annot.protein[annot.protein$symbol == prot,]

# read lod
file_name <- paste0("./QTLscan/addscan_phenotype/", pheno, ".rds")
lod <- readRDS(file=file_name)

# permutation test
perm <- scan1perm(genoprobs = probs,
                  pheno=expr.protein[,prot$id],
                  kinship = K,
                  addcovar = addcovar[,-1],
                  n_perm = 1000,
                  cores = 10,
                  reml = TRUE)

# save perms
file_name <- paste0("./QTLscan/addscan_phenotype/", pheno, "_perm.rds")
saveRDS(perm, file_name)
