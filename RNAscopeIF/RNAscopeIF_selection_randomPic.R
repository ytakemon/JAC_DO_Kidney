# Identify animals that have the average given expression values for the given group.
library(dplyr)
library(ggplot2)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")
load("./shiny_annotation.RData")

# Selection quality:
Gene_name <- "Lrp2"
Data <- "protein"
Sex <- "M"
Age <- c(6, 18)
Num <- 5

samples <- annot.samples

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
# Get expression data
samples$RNAexpr <- expr.mrna[, gene$id]
samples$Protexpr <- expr.protein[, gene$protein_id]

# Subset selction
samples6 <- samples[ samples$Age %in% 6,]
samples18 <- samples[ samples$Age %in% 18,]

samples6F <- samples6[samples6$Sex %in% "F",]
samples6M <- samples6[samples6$Sex %in% "M",]

samples18F <- samples18[samples18$Sex %in% "F",]
samples18M <- samples18[samples18$Sex %in% "M",]

set.seed(1)
sample(samples6F$Mouse.ID, 3)
sample(samples6M$Mouse.ID, 3)
sample(samples18F$Mouse.ID, 3)
sample(samples18M$Mouse.ID, 3)

#> sample(samples6F$Mouse.ID, 3)
#[1] "DO-1252" "DO-1222" "DO-1201"
#> sample(samples6M$Mouse.ID, 3)
#[1] "DO-1041" "DO-1077" "DO-1076"
#> sample(samples18F$Mouse.ID, 3)
#[1] "DO-0781" "DO-0744" "DO-0795"
#> sample(samples18M$Mouse.ID, 3)
#[1] "DO-1007" "DO-1046" "DO-1159"
