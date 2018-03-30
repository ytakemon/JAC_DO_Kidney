#qsub -v script=${R script name} Rsubmit_args.sh
# R/3.4.3
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney.RData")
#library(qtl2)
library(dplyr)

# parameters
addscan.dir <- "./SNPscan/addscansnp_mrna/"
intscan.dir <- "./SNPscan/intscansnp_mrna/"
file.name <- function(i) paste0(annot.mrna$id[i],"_",annot.mrna$symbol[i],".rds")

# define output
output.file <- "./QTLscan/output/scan1snps_eQTLBestperGene.csv"
output <- annot.mrna[,c(1:5,9)]
output$AdditiveSNPid <- output$AdditiveLOD <- output$AdditivePos <-  output$AdditiveChr <- NA
output$FullSNPid <- output$FullLOD <- output$FullPos <-  output$FullChr <- NA
output$AgeIntSNPid <- output$AgeIntLOD <- output$AgeIntPos <-  output$AgeIntChr <- NA

# for
for (i in 1:nrow(output)) {
  #progress
  print("processing", i, "of", nrow(output))

  # load data
  fit_add <- readRDS(paste0(addscan.dir,file.name(i)))
  fit_full <- readRDS(paste0(intscan.dir,file.name(i)))

  # create SNP reference to used for the rest of the analysis
  if(i == 1){
    SNPref <- fit$snpinfo
  }

  # additive
  output$AdditiveLOD[i] <- fit_add$lod[which.max(fit_add$lod)]
  output$AdditiveSNPid[i] <- rownames(fit_add$lod)[which.max(fit_add$lod)]
  output$AdditiveChr[i] <- SNPref %>% filter(snp_id == output$AdditiveSNPid[i]) %>% select(chr)
  output$AdditivePos[i] <- SNPref %>% filter(snp_id == output$AdditiveSNPid[i]) %>% select(pos)

  # Age interactive
  output$FullLOD[i] <- fit_full$lod[which.max(fit_full$lod)]
  output$FullSNPid[i] <- rownames(fit_full$lod)[which.max(fit_full$lod)]
  output$FullChr[i] <- SNPref %>% filter(snp_id == output$FullSNPid[i]) %>% select(chr)
  output$FullPos[i] <- SNPref %>% filter(snp_id == output$FullSNPid[i]) %>% select(pos)

  # Diff
  if(!identical(rownames(fit_add$lod),rownames(fit_full$lod))){
    stop("rsid rownames do not match between fit_add and fit_full!")
  }
  fit_diff <- fit_full$lod - fit_add$lod

  output$AgeIntLOD[i] <- fit_diff[which.max(fit_diff)]
  output$AgeIntSNPid[i] <- rownames(fit_diff)[which.max(fit_diff)]
  output$AgeIntChr[i] <- SNPref %>% filter(snp_id == output$AgeIntSNPid[i]) %>% select(chr)
  output$AgeIntPos[i] <- SNPref %>% filter(snp_id == output$AgeIntSNPid[i]) %>% select(pos)
}
write.csv(output, file=output.file, row.names=FALSE)
