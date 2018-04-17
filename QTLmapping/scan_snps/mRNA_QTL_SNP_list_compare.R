library(tidyverse)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")

select <- 12

eQTL_list <- readr::read_csv("./QTLscan/output/Threshold6_eQTL_intAge.csv", guess_max = 15000) %>% filter(IntAgeChr == select)
SNP_list <- readr::read_csv("./SNPscan/scan1snps_m_diffAgeInt_BestperGene_thr5.csv", guess_max = 4200) %>% filter(AgeIntChr == select)

eQTL_list[eQTL_list$id %in% intersect(eQTL_list$id, SNP_list$id),]



hist(eQTL_list[eQTL_list$id %in% intersect(eQTL_list$id, SNP_list$id),]$IntAgePos * 1e-6)
