#qsub -v script=Paf1_protein_int_permute Rsubmit_args.sh
library(qtl2)
library(qtl2convert)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")
load("./shiny_annotation.RData")
name <- "Paf1"

other.ids <- function(gene.name, level) {
    sel <- which(mRNA.list$symbol == gene.name)[1]
    if (!is.na(sel)) return(mRNA.list[sel,]) else return(c(NA,NA,NA))
  }

id <- other.ids(name)

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "pos")
K <- calc_kinship(probs, "loco", cores=10)
snps$chr <- as.character(snps$chr)
map <- map_df_to_list(map = snps, pos_column = "pos")

addcovar <- model.matrix(~ Sex + Generation + Age + Protein.Batch + Protein.Channel , data = annot.samples)
intcovar <- model.matrix(~ Age, data = annot.samples)

# permutation test
intperm <- scan1perm(genoprobs = probs,
                     pheno=expr.protein[,id$protein_id],
                     kinship = K,
                     addcovar = addcovar[,-1],
                     intcovar = intcovar[,-2],
                     n_perm = 1000,
                     cores = 10,
                     reml = TRUE)

# save perms
file_name <- paste0("./QTLperm/pQTL_Ageint_", name, "_1000perm.rds")
saveRDS(intperm, file_name)
