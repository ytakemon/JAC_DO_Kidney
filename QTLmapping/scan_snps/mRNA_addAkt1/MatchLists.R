library(tidyverse)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")

eQTL_chr12_list <- readr::read_csv("./QTLscan/output/Threshold6_eQTL_intAge.csv", guess_max = 15000) %>% filter(IntAgeChr == "12")
SNP_chr12_list <- readr::read_csv("./SNPscan/scan1snps_m_diffAgeInt_BestperGene_thr5.csv", guess_max = 4200) %>%
  filter(AgeIntChr == "12")

SNP_chr12_list <- SNP_chr12_list %>% mutate(
  eQTL_list = id %in% eQTL_chr12_list$id
)


eQTL_chr12_list$IntAgePos <- eQTL_chr12_list$IntAgePos * 1e-6
