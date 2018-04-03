library(qtl2)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")

not_scanned <- NULL
scanned <- NULL
for (p in 1:nrow(annot.mrna)) {
  print(paste(p, "of ",nrow(annot.mrna)))
  file_name <- paste0("./SNPscan/intscansnp_mrna/", annot.mrna$id[p], "_", annot.mrna$symbol[p], ".rds")
  # in case wall time runs out and the rest need to still be run
  if(file.exists(file_name)){
    scanned <- c(scanned, file_name)
  } else {
    not_scanned <- c(not_scanned, file_name)
  }
}

not_scanned <- NULL
scanned <- NULL
for (p in 1:nrow(annot.protein)) {
  print(paste(p, "of ",nrow(annot.protein)))
  file_name <- paste0("./SNPscan/addscansnp_prot/", annot.protein$id[p], "_", annot.protein$symbol[p], ".rds")
  # in case wall time runs out and the rest need to still be run
  if(file.exists(file_name)){
    scanned <- c(scanned, file_name)
  } else {
    not_scanned <- c(not_scanned, file_name)
  }
}
