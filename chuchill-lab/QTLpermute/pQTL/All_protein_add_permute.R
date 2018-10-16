#qsub -v I=1,script=All_protein_add_permute Rsubmit_args.sh
library(qtl2)
library(qtl2convert)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")
load("./shiny_annotation.RData")
# get input
select <- as.numeric(commandArgs(trailingOnly = TRUE)) # numerical input either 1, 1001, 2001, 3001 ..... last 22,001

# Set range
if ( select == 22000){
  range <- select + 311
} else {
  range <- select + 999
}

# subset for input
subset <- annot.protein[select:range,]

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "pos")
K <- calc_kinship(probs, "loco", cores=10)
snps$chr <- as.character(snps$chr)
map <- map_df_to_list(map = snps, pos_column = "pos")

addcovar <- model.matrix(~ Sex + Generation + Age , data = annot.samples)


for (i in 1:nrow(subset)){
  print(paste("Permuting sequence ", (select - 1) + i))
  id <- subset$id[i]
  symbol <- subset$symbol[i]

  # permutation test
  addperm <- scan1perm(genoprobs = probs,
                       pheno=expr.protein[,id],
                       kinship = K,
                       addcovar = addcovar[,-1],
                       n_perm = 1000,
                       cores = 10,
                       reml = TRUE)

  # save perms
  file_name <- paste0("./QTLperm/pQTL/pQTL_Ageadd_", symbol, "_1000perm.rds")
  saveRDS(addperm, file_name)
}
