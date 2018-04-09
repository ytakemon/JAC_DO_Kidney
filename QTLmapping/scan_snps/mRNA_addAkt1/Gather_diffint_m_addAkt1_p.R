#qsub -v script=${R script name} Rsubmit_args.sh
# R/3.4.3
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney.RData")
options(dplyr.width = Inf)
#library(qtl2)
library(tidyverse)

# parameters
diffscan_dir <- "./SNPscan/diffscansnp_mrna_addAkt1_m"
IntAge_output_file <- "./SNPscan/scan1snps_m_diffAgeInt_addAkt1_m_BestperGene.csv"

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

diffLOD <- gatherdf(dir = diffscan_dir, pattern = "maxLODscan_batch_")

# write table
write_csv(diffLOD, IntAge_output_file)
