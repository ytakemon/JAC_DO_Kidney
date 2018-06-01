library(qtl2)
library(qtl2convert)
library(tidyverse)
#library(ggsci)
options(dplyr.width = Inf) #override column limit
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")
outdir <- "./QTLscan/output/plots/"

# get table to use
IntAge_output_file <- "./QTLscan/scanBestMarker_mrna/BestMarker_BestperGene_mrna.csv"
table <- read_csv(IntAge_output_file)
table12 <- table %>% filter(IntAgeChr == 12, IntAgeLODDiff > 7)
ggplot(table12, aes(x= IntAgePos, y= IntAgeLODDiff)) +
geom_point()

# grab top 4 LODDiff
table12 %>% arrange(-IntAgeLODDiff) %>% select(id, symbol, IntAgeChr, IntAgePos, IntAgeLODDiff)
#   id                 symbol   IntAgeChr IntAgePos IntAgeLODDiff
#   <chr>              <chr>    <chr>         <dbl>         <dbl>
# 1 ENSMUSG00000036850 Mrpl41   12             87.1          11.5
# 2 ENSMUSG00000001833 Sept7    12             89.0          11.0
# 3 ENSMUSG00000042380 Smim12   12            103.           11.0
# 4 ENSMUSG00000090306 Adh6-ps1 12            103.           10.6
# 5 ENSMUSG00000034118 Tpst1    12             89.6          10.6

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

Sept7 <- MakeEachQTL("ENSMUSG00000001833", "Sept7", "eQTL", c(0,10))
#show plots
pdf(paste0(outdir,"Sept7_eQTL_6mo.pdf"), height = 4, width = 8)
Sept7[[1]]
dev.off()

pdf(paste0(outdir,"Sept7_eQTL_12mo.pdf"), height = 4, width = 8)
Sept7[[2]]
dev.off()

pdf(paste0(outdir,"Sept7_eQTL_18mo.pdf"), height = 4, width = 8)
Sept7[[3]]
dev.off()

Smim12 <- MakeEachQTL("ENSMUSG00000042380", "Smim12", "eQTL", c(0,10))
#show plots
pdf(paste0(outdir,"Smim12_eQTL_6mo.pdf"), height = 4, width = 8)
Smim12[[1]]
dev.off()

pdf(paste0(outdir,"Smim12_eQTL_12mo.pdf"), height = 4, width = 8)
Smim12[[2]]
dev.off()

pdf(paste0(outdir,"Smim12_eQTL_18mo.pdf"), height = 4, width = 8)
Smim12[[3]]
dev.off()

Adh6 <- MakeEachQTL("ENSMUSG00000090306", "Adh6-ps1", "eQTL", c(0,8))
#show plots
pdf(paste0(outdir,"Adh6-ps1_eQTL_6mo.pdf"), height = 4, width = 8)
Adh6[[1]]
dev.off()

pdf(paste0(outdir,"Adh6-ps1_eQTL_12mo.pdf"), height = 4, width = 8)
Adh6[[2]]
dev.off()

pdf(paste0(outdir,"Adh6-ps1_eQTL_18mo.pdf"), height = 4, width = 8)
Adh6[[3]]
dev.off()

Tpst1 <- MakeEachQTL("ENSMUSG00000034118", "Tpst1", "eQTL", c(0,12))
#show plots
pdf(paste0(outdir,"Tpst1_eQTL_6mo.pdf"), height = 4, width = 8)
Tpst1[[1]]
dev.off()

pdf(paste0(outdir,"Tpst1_eQTL_12mo.pdf"), height = 4, width = 8)
Tpst1[[2]]
dev.off()

pdf(paste0(outdir,"Tpst1_eQTL_18mo.pdf"), height = 4, width = 8)
Tpst1[[3]]
dev.off()

#   id                 symbol   IntAgeChr IntAgePos IntAgeLODDiff
#   <chr>              <chr>    <chr>         <dbl>         <dbl>
# 1 ENSMUSG00000036850 Mrpl41   12             87.1          11.5
# 2 ENSMUSG00000001833 Sept7    12             89.0          11.0
# 3 ENSMUSG00000042380 Smim12   12            103.           11.0
# 4 ENSMUSG00000090306 Adh6-ps1 12            103.           10.6
# 5 ENSMUSG00000034118 Tpst1    12             89.6          10.6
