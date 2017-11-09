setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")

akt <- read.csv("./QTLscan/output/eQTLintAkt1thr6_chr12.csv", stringsAsFactors = FALSE)
AKT <- read.csv("./QTLscan/output/eQTLintAKT1thr6_chr12.csv", stringsAsFactors = FALSE)
AKT_wb_tot <- read.csv("./QTLscan/output/eQTLintPanTotAkt1thr6_chr12.csv", stringsAsFactors = FALSE)

# Identify gene with a LOD drop of 2
LODdrop <- 2

akt$drop <- akt$IntAgeLODDiff - akt$addIntAgeLODDiff
AKT$drop <- AKT$IntAgeLODDiff - AKT$addIntAgeLODDiff
AKT_wb_tot$drop <- AKT_wb_tot$IntAgeLODDiff - AKT_wb_tot$addIntAgeLODDiff

akt$drop2 <- akt$drop > 2
AKT$drop2 <- AKT$drop > 2
AKT_wb_tot$drop2 <- AKT_wb_tot$drop > 2

akt <- akt[akt$drop2 == TRUE,]
AKT <- AKT[AKT$drop2 == TRUE,]
AKT_wb_tot <- AKT_wb_tot[AKT_wb_tot$drop2 == TRUE,]

> nrow(akt)
[1] 186  11
> nrow(AKT)
[1] 14 11
> nrow(AKT_wb_tot)
[1] 209  11

# mRNA - shotgun proteomics
if (nrow(akt) > nrow(AKT)){
  mRNA_prot <- akt[akt$id %in% AKT$id, ]
  mRNA_prot$pair <- "mRNA_prot"
} else {
  mRNA_prot <- AKT[AKT$id %in% akt$id, ]
  mRNA_prot$pair <- "mRNA_prot"
}
# mRNA - WB Pan total
if (nrow(akt) > nrow(AKT_wb_tot)){
  mRNA_WBPanTot <- akt[akt$id %in% AKT_wb_tot$id,]
  mRNA_WBPanTot$pair <- "mRNA_WBPanTot"
} else {
  mRNA_WBPanTot <- AKT_wb_tot[AKT_wb_tot$id %in% akt$id,]
  mRNA_WBPanTot$pair <- "mRNA_WBPanTot"
}
# shotgun proteomics - WB Pan total
if (nrow(AKT) > nrow(AKT_wb_tot)){
  prot_WBPanTot <- AKT[AKT$id %in% AKT_wb_tot$id,]
  prot_WBPanTot$pair <- "prot_WBPanTot"
} else {
  prot_WBPanTot <- AKT_wb_tot[AKT_wb_tot$id %in% AKT$id,]
  prot_WBPanTot$pair <- "prot_WBPanTot"
}

combine <- rbind(mRNA_prot, mRNA_WBPanTot, prot_WBPanTot)
write.csv(combine, file = "./QTLscan/output/Mediation_common_genes.csv", row.names = FALSE)

# > nrow(mRNA_prot)
# [1] 5
# > nrow(mRNA_WBPanTot)
# [1] 209
# > nrow(prot_WBPanTot)
# [1] 5
