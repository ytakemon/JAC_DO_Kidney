plist <- as.numeric(commandArgs(trailingOnly = TRUE))

cat(paste("Mapping", length(plist), "genes", "\n"))
print(plist)

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

plist <- plist[plist<=ncol(expr.mrna)]

for (p in plist) {

  cat("Scanning ",which(p==plist)," out of ",length(plist),"\n")

  addcovar <- model.matrix(~ Sex + Age + Generation, data=annot.samples)

  # lod score
  lod <- scan1snps(genoprobs=probs,
               kinship=Glist,
               pheno=expr.mrna[,p],
               addcovar=addcovar[,-1], cores=10, reml=TRUE)

  # save lod object
  file_name <- paste0("./QTLscan/addscan_mrna/", annot.mrna$id[p], "_", annot.mrna$symbol[p], ".rds")
  saveRDS(lod, file=file_name)
}

# load example data and calculate genotype probabilities
     file <- paste0("https://raw.githubusercontent.com/rqtl/",
                    "qtl2data/master/DOex/DOex.zip")
     DOex <- read_cross2(file)
     probs <- calc_genoprob(DOex, error_prob=0.002)

     snpdb_file <- system.file("extdata", "cc_variants_small.sqlite", package="qtl2")
     queryf <- create_variant_query_func(snpdb_file)

     out <- scan1snps(probs, DOex$pmap, DOex$pheno, query_func=queryf, chr=2, start=97, end=98)
     ## End(Not run)
