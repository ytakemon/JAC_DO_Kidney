library(qtl2)
library(qtl2convert)
library(tidyverse)
#library(ggsci)
options(dplyr.width = Inf) #override column limit
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")
outdir <- "./QTLscan/output/plots/"

snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="20"] <- "X"
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
map <- map_df_to_list(map = snps, pos_column = "bp")
K <- calc_kinship(probs, type = "loco", cores = 10)

# visualize
MakeEachQTL <- function(EnsID, symbol, type, ylim){

  # remove after function validated
  #EnsID <- pick$id
  #symbol <- pick$symbol
  #type <- "eQTL"

  # identify type and define path
  if (type == "eQTL"){
    path <- "./QTLscan/addscan_mrna_perAge"
  } else if (type == "pQTL"){
    path <- "./QTLscan/addscan_prot_perAge"
  } else {
    stop("don't recognise type, use either eQTL or pQTL")
  }

  # get scan1 from each time point
  qtl_6mo <- readRDS(paste0(path,"/6mo/",EnsID,"_",symbol,".rds"))
  qtl_12mo <- readRDS(paste0(path,"/12mo/",EnsID,"_",symbol,".rds"))
  qtl_18mo <- readRDS(paste0(path,"/18mo/",EnsID,"_",symbol,".rds"))

  #plotting (not sure how to do this quitely in the background)
  plot_scan1(qtl_6mo, map, ylim = ylim)
  title(paste(symbol, "(", EnsID, ")", "@6mo"))
  p6 <- recordPlot()
  plot_scan1(qtl_12mo, map, ylim = ylim)
  title(paste(symbol, "(", EnsID, ")", "@12mo"))
  p12 <- recordPlot()
  plot_scan1(qtl_18mo, map, ylim = ylim)
  title(paste(symbol, "(", EnsID, ")", "@18mo"))
  p18 <- recordPlot()

  allplots <- list(p6,p12,p18)
  return(allplots)
  # see using:
  # allplots[[1]] for 6 months
  # allplots[[2]] for 12 months
  # allplots[[3]] for 18 months
}

Mrpl <- MakeEachQTL("ENSMUSG00000036850", "Mrpl41", "eQTL", c(0,11.5))
#show plots
pdf(paste0(outdir,"Mrpl_eQTL_6mo.pdf"), height = 4, width = 8)
Mrpl[[1]]
dev.off()

pdf(paste0(outdir,"Mrpl_eQTL_12mo.pdf"), height = 4, width = 8)
Mrpl[[2]]
dev.off()

pdf(paste0(outdir,"Mrpl_eQTL_18mo.pdf"), height = 4, width = 8)
Mrpl[[3]]
dev.off()
