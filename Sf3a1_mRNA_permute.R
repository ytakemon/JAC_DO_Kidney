#qsub -v script=Sf3a1_mRNA_permute Rsubmit_args.sh
library(qtl2)
library(qtl2convert)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")
load("./shiny_annotation.RData")
name <- "Sf3a1"

other.ids <- function(gene.name, level) {
    sel <- which(mRNA.list$symbol == gene.name)[1]
    if (!is.na(sel)) return(mRNA.list[sel,]) else return(c(NA,NA,NA))
  }

id <- other.ids(name)

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "pos")
K <- calc_kinship(probs, "loco", cores=3)
snps$chr <- as.character(snps$chr)
map <- map_df_to_list(map = snps, pos_column = "pos")

addcovar <- model.matrix(~ Sex + Generation + Age , data = annot.samples)

# permutation test
addperm <- scan1perm(genoprobs = probs,
                     pheno=expr.mrna[,id$id],
                     kinship = K,
                     addcovar = addcovar[,-1],
                     n_perm = 1000,
                     cores = 10,
                     reml = TRUE)

# save perms
file_name <- paste0("./QTLperm/eQTL_", name, "_1000perm.rds")
saveRDS(perm, file_name)
