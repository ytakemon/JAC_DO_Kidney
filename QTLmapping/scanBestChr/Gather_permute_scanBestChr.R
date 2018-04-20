setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
#load("./RNAseq_data/DO188b_kidney.RData")
options(dplyr.width = Inf)
#library(qtl2)
library(tidyverse)

# parameters
# mrna
diradd_m <- "./QTLscan/scanBestChr_mrna/perm/add"
dirfull_m <- "./QTLscan/scanBestChr_mrna/perm/full"
outputadd_m <- "./QTLscan/scanBestChr_mrna/perm/All_addperm_m.csv"
outputfull_m <- "./QTLscan/scanBestChr_mrna/perm/All_Fullperm_m.csv"
# protein
diradd_p <- "./QTLscan/scanBestChr_protein/perm/add"
dirfull_p <- "./QTLscan/scanBestChr_protein/perm/full"
outputadd_p <- "./QTLscan/scanBestChr_protein/perm/All_addperm_p.csv"
outputfull_p <- "./QTLscan/scanBestChr_protein/perm/All_Fullperm_p.csv"

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

# mrna
Permadd_m <- gatherdf(dir = diradd_m, pattern = "perm10scan_batch")
Permfull_m <- gatherdf(dir = dirfull_m, pattern = "perm10scan_batch")
# protein
Permadd_p <- gatherdf(dir = diradd_p, pattern = "perm10scan_batch")
Permfull_p <- gatherdf(dir = dirfull_p, pattern = "perm10scan_batch")

process_perm <- function(perm){
  perm <- select(perm,contains("perm"))
  print(quantile(as.matrix(perm),c(0.6,0.8,0.9,0.95),na.rm = TRUE))
  hist(as.matrix(perm),100)
  abline(v = quantile(as.matrix(perm),c(0.6,0.8,0.9,0.95),na.rm=TRUE))
}

process_perm(Permadd_m)
process_perm(Permfull_m)
process_perm(Permadd_p)
process_perm(Permfull_p)

# write table
write_csv(Permadd_m, outputadd_m)
write_csv(Permfull_m, outputfull_m)
write_csv(Permadd_p, outputadd_p)
write_csv(Permfull_p, outputfull_p)
