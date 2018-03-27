plist <- as.numeric(commandArgs(trailingOnly = TRUE))

cat(paste("Mapping", length(plist), "genes", "\n"))
print(plist)

library(qtl2)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")
snpdb <- "/hpcdata/gac/resource/CCsnps/cc_variants_v2.sqlite"

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
map <- map_df_to_list(map = snps, pos_column = "bp")
query_func <- create_variant_query_func(snpdb)

plist <- plist[plist<=ncol(expr.mrna)]

for (p in plist) {

  cat("Scanning ",which(p==plist)," out of ",length(plist),"\n")

  addcovar <- model.matrix(~ Sex + Age + Generation, data=annot.samples)

  # Perform scan1snps
  # If query_func is given, but start and end are empty, it should calcualte for all chromosomes.
  snpsOut <- scan1snps(genoprobs=probs,
               map = map,
               kinship=Glist,
               pheno=expr.mrna[,p],
               addcovar=addcovar[,-1],
               query_func = query_func,
               chr =c(1:19,"X"),
               start = 0,
               end = 200,
               keep_all_snps = TRUE,
               cores=20, reml=TRUE)

  # save lod object
  file_name <- paste0("./SNPscan/addscansnp_mrna/", annot.mrna$id[p], "_", annot.mrna$symbol[p], ".rds")
  saveRDS(snpsOut, file=file_name)
}
