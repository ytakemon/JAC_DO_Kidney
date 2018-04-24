setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney.RData")
options(dplyr.width = Inf)
#library(qtl2)
library(tidyverse)
library(qtl2)
library(qtl2convert)

# read in data
# age interactive gene list (mrna)
infile <- "./QTLscan/scanBestMarker_mrna/BestMarker_BestperGene_mrna.csv"
data <- readr::read_csv(infile)
data_thr8_m <- filter(data, IntAgeChr == "12", IntAgeLODDiff > 8, ) %>% arrange(-IntAgeLODDiff)

#write_csv(data_thr8_m, "./QTLscan/scanBestMarker_mrna/BestMarker_BestperGene_mrna_thr8.csv")

# age interactive gene list (protein)
infile <- "./QTLscan/scanBestMarker_protein/BestMarker_BestperGene_protein.csv"
data <- readr::read_csv(infile)
data_thr10_p <- filter(data, IntAgeChr == "7", IntAgeLODDiff > 10) %>% arrange(-IntAgeLODDiff)
#write_csv(data_thr10_p, "./QTLscan/scanBestMarker_protein/BestMarker_BestperGene_protein_thr10.csv")

# prepare for qtl plot
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")

# plotting function
QTLplot <- function(id,symbol,type,model){
  if(type == "mrna"){
    out <- "./QTLscan/scanBestMarker_mrna/QTLplotValid/"
    if(model == "add"){  # The data are the same as the by chr scans
      dir <- "./QTLscan/scanBestChr_mrna/addscan/"
    } else if(model == "full"){
      dir <- "./QTLscan/scanBestChr_mrna/fullscan/"
    }
  } else if (type == "protein"){
    out <- "./QTLscan/scanBestMarker_protein/QTLplotValid/"
    if(model == "add"){ # The data are the same as the by chr scans
      dir <- "./QTLscan/scanBestChr_protein/addscan/"
    } else if (model == "full"){
      dir <- "./QTLscan/scanBestChr_protein/fullscan/"
    }
  }

  #find and get
  file <- list.files(path = dir, pattern = paste0(id,"_",symbol,".rds"), full.names = TRUE)
  scan <- readRDS(file)

  #plot
  pdf(paste0(out,"QTLplot_",symbol,"_",type,"_",model,".pdf"),height=3,width=7)
  plot_scan1(scan, map, ylim = c(0,16))
  title(paste0(id," ",symbol, " (",type,") ", model,"scan model" ))
  dev.off()
}

QTLplot(data_thr8_m$id[1],data_thr8_m$symbol[1], "mrna", "add")
QTLplot(data_thr8_m$id[1],data_thr8_m$symbol[1], "mrna", "full")
QTLplot(data_thr8_m$id[2],data_thr8_m$symbol[2], "mrna", "add")
QTLplot(data_thr8_m$id[2],data_thr8_m$symbol[2], "mrna", "full")
QTLplot(data_thr8_m$id[3],data_thr8_m$symbol[3], "mrna", "add")
QTLplot(data_thr8_m$id[3],data_thr8_m$symbol[3], "mrna", "full")
QTLplot(data_thr8_m$id[4],data_thr8_m$symbol[4], "mrna", "add")
QTLplot(data_thr8_m$id[4],data_thr8_m$symbol[4], "mrna", "full")

QTLplot(data_thr10_p$id[1],data_thr10_p$symbol[1], "protein", "add")
QTLplot(data_thr10_p$id[1],data_thr10_p$symbol[1], "protein", "full")
QTLplot(data_thr10_p$id[2],data_thr10_p$symbol[2], "protein", "add")
QTLplot(data_thr10_p$id[2],data_thr10_p$symbol[2], "protein", "full")
QTLplot(data_thr10_p$id[3],data_thr10_p$symbol[3], "protein", "add")
QTLplot(data_thr10_p$id[3],data_thr10_p$symbol[3], "protein", "full")
QTLplot(data_thr10_p$id[4],data_thr10_p$symbol[4], "protein", "add")
QTLplot(data_thr10_p$id[4],data_thr10_p$symbol[4], "protein", "full")
