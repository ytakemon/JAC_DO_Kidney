library(ggplot2)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")
load("./shiny_annotation.RData")

# Selection quality:
Gene_name <- "Lrp2"
Data <- "mRNA"
Sex <- "M"
Age <- c(6, 18)
Num <- 3 # Number of animals to select

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
gene <- other.ids(Gene_name, Data)

if(Data == "mRNA"){
  annot.samples$expr <- expr.mrna[, gene$id]
} else if (Data == "protein"){
  annot.samples$expr <- expr.protein[, gene$protein_id]
}
