# run with:
# qsub -v script=pQTL_prot_pERK1Upstream_int Rsubmit_args.sh
library(qtl2)
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO188b_kidney.RData")

# Get query list and see if they exist in RNA-seq or Shotgun-proteomics
query_list <- read.csv("./SpecificQ/pERK1_addtotERK1_LODpeak_genes.csv", header = TRUE, stringsAsFactors = FALSE)
E_query_list <- query_list[query_list$ensembl_gene_id %in% annot.mrna$id,]
P_query_list <- query_list[query_list$ensembl_gene_id %in% annot.protein$gene_id,]

# Get list of genes with trans pQTL
list <- read.csv("./QTLscan/output/Threshold6_pQTL_intAge_pbatch.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == 7, ]
list <- list$id # ENS protein id

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "pos")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "pos")

# Get mediator from query_list[i] to annotation
for (g in 1:nrow(P_query_list)){

  med_query <- P_query_list$ensembl_gene_id[g]
  med_query <- annot.protein[annot.protein$gene_id == med_query,]$id
  annot.samples$med_query <- expr.protein[,med_query]

  for (p in 1:length(list)){

    cat("Scanning ", p, " out of ", length(list), ", of query add list ", g, " out of ", nrow(P_query_list), "\n")

    addcovar <- model.matrix(~ Sex + Age + Generation + Protein.Batch + Protein.Channel + med_query,
                            data=annot.samples)
    intcovar <- model.matrix(~ Age, data=annot.samples)

    p <- list[p]
    p <- annot.protein[annot.protein$id == p, ]

    # Lod score
    lod <- scan1(genoprobs=probs,
                 kinship=Glist,
                 pheno=expr.protein[, p$id],
                 addcovar=addcovar[,-1],
                 intcovar=intcovar[,-1],
                 cores=10, reml=TRUE)

    # Save lod object
    file_name <- paste0("./QTLscan/intscan_prot_pERK1Upstream/Add_", med_query,"_Medscan_", p$id, "_", p$symbol, ".rds")
    saveRDS(lod, file = file_name)

  }
}
