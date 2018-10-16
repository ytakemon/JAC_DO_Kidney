#qsub -v script=${R script name} Rsubmit_args.sh
# R/3.4.3
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney.RData")
options(dplyr.width = Inf)
#library(qtl2)
library(tidyverse)

# parameters
diffscan_dir <- "./QTLscan/scanBestMarker_SexInt_protein/maxLODscan"
IntSex_output_file <- "./QTLscan/scanBestMarker_SexInt_protein/BestMarker_SexInt_BestperGene_protein.csv"

gatherdf <- function(dir, pattern){
  file_list <- list.files(dir,pattern,full.names=TRUE)

  #read and gather
  output <- NULL
  for(i in file_list){
    one <- readr::read_csv(i) %>% mutate(chr = as.character(chr))
    output <- bind_rows(output,one)
  }

  # return df
  output <- output %>% arrange(id)
  return(output)
}

diffLOD <- gatherdf(dir = diffscan_dir, pattern = "maxLODscan_batch_")

# Oops I have to rename the Interaction to Sex!
diffLOD <- diffLOD %>% rename(IntSexChr = IntAgeChr,
  IntSexPos = IntAgePos,
  IntSexLODDiff = IntAgeLODDiff)

# write table
write_csv(diffLOD, IntSex_output_file)
