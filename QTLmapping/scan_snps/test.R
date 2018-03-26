# qsub -v script=test Rsubmit_args.sh
library(qtl2)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")
snpdb <- "/hpcdata/gac/resource/CCsnps/cc_variants_v2.sqlite"

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "pos")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "pos")
query_func <- create_variant_query_func(snpdb)

p <- 1

  addcovar <- model.matrix(~ Sex + Age + Generation, data=annot.samples)

  # Perform scan1snps
  # If query_func is given, but start and end are empty, it should calcualte for all chromosomes.
  snpsOut <- scan1snps(genoprobs=probs,
               map = map,
               kinship=Glist,
               pheno=expr.mrna[,p],
               addcovar=addcovar[,-1],
               query_func = query_func,
               keep_all_snps = TRUE,
               cores=20, reml=TRUE)

  # save lod object
  file_name <- paste0("./SNPscan/addscansnp_mrna/", annot.mrna$id[p], "_", annot.mrna$symbol[p], ".rds")
  saveRDS(snpsOut, file=file_name)
