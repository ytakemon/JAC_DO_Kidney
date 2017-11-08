# run with:
# qsub -v script=eQTL_Akt1Add_WbPanTotAkt2 Rsubmit_args.sh
# qsub -v script=eQTL_Akt1Add_WbPanTotAkt4 Rsubmit_args.sh
# qsub -v script=eQTL_Akt1Add_WbPanTotAkt6 Rsubmit_args.sh
# qsub -v script=eQTL_Akt1Add_WbPanTotAkt8 Rsubmit_args.sh
# qsub -v script=eQTL_Akt1Add_WbPanTotAkt9 Rsubmit_args.sh

library(qtl2)
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO188b_kidney.RData")
akt <- read.delim("./Phenotype/phenotypes/JAC_WB_kidney_AKT.txt")

# Note from George: Multiple bands in Akt are not know to have different
# functions, so Pan-AKT was created to reflect total expression.
# Will use Pan-AKT data for this analysis.
akt <- akt[complete.cases(akt$Total_Pan.AKT),]

# cleanup Akt data
akt$duplicated <- (duplicated(akt$ID) | duplicated(akt$ID, fromLast = TRUE))
akt <- akt[akt$duplicated == FALSE,]
akt <- akt[akt$ID %in% rownames(genoprobs),] # 159 animals
rownames(akt) <- akt$ID

# Get list of genes with trans eQTL
list <- read.csv("./QTLscan/output/Threshold6_eQTL_intAge.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == 12, ]
list <- list$symbol

# match animal data
genoprobs <- genoprobs[rownames(genoprobs) %in% rownames(akt),,]
annot.samples <- annot.samples[rownames(annot.samples) %in% rownames(akt),]
expr.mrna <- expr.mrna[rownames(expr.mrna) %in% rownames(akt),]

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")
K <- calc_kinship(probs, type = "loco")

# Identify gene name
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

# this piece should can be broken apart if the list is long...
# doesn't apply for interactive QLT run.
# eg.
#for (p in 1:200){
#for (p in 201:400){
#for (p in 401:600){
#for (p in 601:800){
#for (p in 801:length(list)){

for (p in 1:length(list)){

  cat("Scanning ",p," out of ",length(list),"\n")
  addcovar <- model.matrix(~ Sex + Age + Generation + log(akt$Total_Pan.AKT), data=annot.samples)

  p <- list[p]
  p <- other.ids(p, "mRNA")

  # Lod score
  lod <- scan1(genoprobs=probs,
               kinship=K,
               pheno=expr.mrna[, p$id],
               addcovar=addcovar[,-1],
               cores=10, reml=TRUE)

  # Save lod object
  file_name <- paste0("./QTLscan/addscan_mrna_WB_PanTotAkt1/", p$id, "_", p$symbol, ".rds")
  saveRDS(lod, file = file_name)
}
