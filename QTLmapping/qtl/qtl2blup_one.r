number_of_parts <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
which_part <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

library(qtl2geno)
library(qtl2scan)
library(qtl2convert)

# load kidney data
load("DO192_kidney_gary.RData")

# define covariates 
covar <- model.matrix(~ Age*Sex + Generation, data=pheno)[,-1]

# covert genotypes to qtl2 format
genoprobs <- probs_doqtl_to_qtl2(probs, snps, pos_column = "bp")
map <- map_df_to_list(map = snps, pos_column = "bp")

# generate table with all scans
table_of_scans <- rbind(data.frame(pheno="expr.mrna", col = 1:ncol(expr.mrna), addcovar="covar", 
                                   genoprobs="genoprobs", kinship = "Glist", stringsAsFactors = FALSE),
                        data.frame(pheno="expr.protein", col = 1:ncol(expr.protein), addcovar="covar", 
                                   genoprobs="genoprobs", kinship = "Glist", stringsAsFactors = FALSE))

# split N into k subsets and report j-th range
split_set <- function(N,k,j) {
  Nsubset <- ceiling(N / k)
  j_th_from <- Nsubset * (j-1) + 1
  j_th_to <- min(Nsubset * j, N)
  return(c(j_th_from, j_th_to))
}

from_i = split_set(nrow(table_of_scans), number_of_parts, which_part)[1]
to_i = split_set(nrow(table_of_scans), number_of_parts, which_part)[2]
lods <- list()
effects <- list()
chrs <- c(as.character(1:19), "X")

for (i in from_i:to_i) {

  cat("Scanning ",i-from_i+1," out of ",to_i-from_i+1,"\n")

  tmp <- NULL
  tmp <- scan1(genoprobs=get(table_of_scans$genoprobs[i]),
               kinship=get(table_of_scans$kinship[i]),
               pheno=get(table_of_scans$pheno[i])[, table_of_scans$col[i]],
               addcovar=get(table_of_scans$addcovar[i]), cores=1)
  lods[[i]] <- tmp$lod
  
  tmp_coefs <- NULL
  for (c in chrs) {
    tmp2 <- scan1blup(genoprobs = get(table_of_scans$genoprobs[i])[,c], 
                      pheno = get(table_of_scans$pheno[i])[, table_of_scans$col[i]], 
                      addcovar = get(table_of_scans$addcovar[i]), 
                      kinship=get(table_of_scans$kinship[i])[[c]])
    tmp_coefs <- rbind(tmp_coefs, tmp2$coef[,LETTERS[1:8]])
  }  
  effects[[i]] <- tmp_coefs
}  

file.name1 <- paste0("qtl2blup/lods_", number_of_parts, "_", which_part, ".rds")
saveRDS(lods, file=file.name1)
file.name2 <- paste0("qtl2blup/effects_", number_of_parts, "_", which_part, ".rds")
saveRDS(effects, file=file.name2)

