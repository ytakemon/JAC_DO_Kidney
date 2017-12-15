# run with:
# qsub -v script=pQTL_ERK1PhosAdd Rsubmit_args.sh
library(qtl2)
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO188b_kidney.RData")

# Get query list and see if they exist in RNA-seq or Shotgun-proteomics
query_list <- read.csv("./SpecificQ/pERK1_addtotERK1_LODpeak_genes.csv")
E_query_list <- query_list[query_list$ensembl_gene_id %in% annot.mrna$id,]
P_query_list <- query_list[query_list$ensembl_gene_id %in% annot.protein$gene_id,]

# Get list of genes with trans pQTL
list <- read.csv("./QTLscan/output/Threshold6_pQTL_intAge_pbatch.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == 7, ]
list <- list$id # ENS protein id

# add query_list[i] to annotation
for (g in 1:)
annot.samples$erk1 <- erk1[,pheno]

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")

for (p in 1:length(list)){

  cat("Scanning ",p," out of ",length(list),"\n")
  P_addcovar <- model.matrix(~ Sex + Age + Generation + Protein.Batch + Protein.Channel + log(erk1), data=annot.samples)

  p <- list[p]
  p <- annot.protein[annot.protein$id == p,]

  # Lod score
  lod <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=expr.protein[, p$id],
               addcovar=addcovar[,-1],
               cores=10, reml=TRUE)

  # Save lod object
  file_name <- paste0("./QTLscan/addscan_pERK1_upstream/", p$id, "_", p$symbol, ".rds")
  saveRDS(lod, file = file_name)
}
