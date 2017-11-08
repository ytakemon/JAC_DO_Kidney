# run with:
# qsub -v script=eQTL_Akt1Add2 Rsubmit_args.sh
# qsub -v script=eQTL_Akt1Add4 Rsubmit_args.sh
# qsub -v script=eQTL_Akt1Add6 Rsubmit_args.sh
# qsub -v script=eQTL_Akt1Add8 Rsubmit_args.sh
# qsub -v script=eQTL_Akt1Add9 Rsubmit_args.sh

library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO188b_kidney.RData")

# Get list of genes with trans eQTL
list <- read.csv("./QTLscan/output/Threshold6_eQTL_intAge.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == 12, ]
list <- list$symbol

# Identify gene name
trans <- "Akt1"
other.ids <- function(gene.name, level) {
  if (level == "mRNA") {
    sel <- which(mRNA.list$symbol == gene.name)[1]
    if (!is.na(sel)) return(mRNA.list[sel,]) else return(c(NA,NA,NA))
  }
  if (level == "protein") {
    sel <- which(protein.list$symbol == gene.name)[1]
    if (!is.na(sel)) return(protein.list[sel,]) else return(c(NA,NA,NA))
  }
}
trans <- other.ids(trans, "mRNA")

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")

# this piece should can be broken apart if the list is long...
# doesn't apply for interactive QLT run.
# eg.
#for (p in 1:200){
#for (p in 201:400){
#for (p in 401:600){
#for (p in 601:800){
#for (p in 801:length(list)){=

for (p in 1:length(list)){

  cat("Scanning ",p," out of ",length(list),"\n")
  addcovar <- model.matrix(~ Sex + Age + Generation + expr.protein[,trans$protein_id], data=annot.samples)

  p <- list[p]
  p <- other.ids(p, "mRNA")

  # Lod score
  lod <- scan1(genoprobs=probs,
               kinship=Glist,
               pheno=expr.mrna[, p$id],
               addcovar=addcovar[,-1],
               cores=10, reml=TRUE)

  # Save lod object
  file_name <- paste0("./QTLscan/addscan_mrna_AKT1/", p$id, "_", p$symbol, ".rds")
  saveRDS(lod, file = file_name)
}
