#qsub -v script=${R script name} Rsubmit_args.sh
# R/3.4.3
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney.RData")
options(dplyr.width = Inf)
#library(qtl2)
library(tidyverse)

# parameters
addscan_dir <- "./SNPscan/addscansnp_prot"
intscan_dir <- "./SNPscan/intscansnp_prot"
Add_output_file <- "./SNPscan/scan1snps_p_Add_BestperGene.csv"
IntAge_output_file <- "./SNPscan/scan1snps_p_AgeInt_BestperGene.csv"

gatherdf <- function(dir, pattern){
  file_list <- list.files(dir,pattern,full.names=TRUE)

  #read and gather
  output <- NULL
  for(i in file_list){
    one <- read_csv(i) %>% mutate(chr = as.character(chr))
    output <- bind_rows(output,one)
  }

  # return df
  output <- output %>% arrange(id)
  return(output)
}

addLOD <- gatherdf(dir = addscan_dir, pattern = "maxLODscan_batch_")
fullLOD <- gatherdf(dir = intscan_dir, pattern = "maxLODscan_batch_")

# combine
output <- addLOD %>% mutate(
  FullChr = fullLOD$FullChr,
  FullPos = fullLOD$FullPos,
  FullLOD = fullLOD$FullLOD,
  IntAgeChr = AdditiveChr,
  IntAgePos = AdditivePos,
  IntAgeLODDiff = FullLOD - AdditiveLOD
)

# get relevant columns
output_add <- output %>% select(id, gene_id, symbol, chr, start, end, strand, biotype,
  AdditiveChr, AdditivePos, AdditiveLOD)
output_ageint <- output %>% select(id, gene_id, symbol, chr, start, end, strand, biotype,
  IntAgeChr, IntAgePos, FullLOD, IntAgeLODDiff)

# write table
write_csv(output_add, Add_output_file)
write_csv(output_ageint, IntAge_output_file)
