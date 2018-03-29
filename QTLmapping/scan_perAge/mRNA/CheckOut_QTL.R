library(qtl2)
library(qtl2convert)
library(tidyverse)
library(grid)
#library(ggsci)
options(dplyr.width = Inf) #override column limit
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")

# load eQTL data
eQTL_best_6 <- read_csv("./QTLscan/output/6mo_eQTLBestperGene.csv")
eQTL_best_12 <- read_csv("./QTLscan/output/12mo_eQTLBestperGene.csv")
eQTL_best_18 <- read_csv("./QTLscan/output/18mo_eQTLBestperGene.csv")
# load pQTL data
pQTL_best_6 <- read_csv("./QTLscan/output/6mo_pQTLBestperGene.csv")
pQTL_best_12 <- read_csv("./QTLscan/output/12mo_pQTLBestperGene.csv")
pQTL_best_18 <- read_csv("./QTLscan/output/18mo_pQTLBestperGene.csv")

# create function to generate df that made eQTL plots (incase)
# no threshold is set at the moment
Getgmb <- function(eQTL_best){
  chr_full <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X","Y", "MT")
  AddQTL <- eQTL_best %>%
              mutate(
                chr = factor(chr, levels = chr_full),
                AdditiveChr = factor(AdditiveChr, levels= chr_full)
              ) %>%
              filter(!(chr %in% c("Y","MT")))

  # Annotate postion with genes and save file for sharing
  save <- arrange(AddQTL, AdditiveChr, AdditivePos)
  # save annotated list for sharing
  #write.csv(save, "./QTLscan/output/Threshold6_eQTL_6mo.csv", row.names = FALSE, quote = FALSE)


  # Convert transcript and qtl position relative to chromosome positions
  # Convert to megabases
  chrlen <- sapply(split(AddQTL$end, AddQTL$chr), max) * 1e-6
  chrsum <- c(0, cumsum(chrlen))
  names(chrsum) = names(chrlen)

  t.gmb <- AddQTL$start * 1e-6 # Transcript
  q.gmb <- AddQTL$AdditivePos * 1e-6 # qtl

  # Cumulative sum of previosu positions
  for(i in c(2:19, "X")) {

    wh <- which(AddQTL$chr == i)
    t.gmb[wh] <- t.gmb[wh] + chrsum[i]

    wh <- which(AddQTL$AdditiveChr == i)
      q.gmb[wh] <- q.gmb[wh] + chrsum[i]
  }
  AddQTL$t_gbm <- t.gmb
  AddQTL$q_gbm <- q.gmb
  return(AddQTL)
}

# Run function of each age
# eQTL
eQTL_6mo <- Getgmb(eQTL_best_6)
eQTL_12mo <- Getgmb(eQTL_best_12)
eQTL_18mo <- Getgmb(eQTL_best_18)

# pQTL
pQTL_6mo <- Getgmb(pQTL_best_6)
pQTL_12mo <- Getgmb(pQTL_best_12)
pQTL_18mo <- Getgmb(pQTL_best_18)

# Compare two time points
CompQTL <- function(qtl1, qtl2){

  # remove from function once done
  #qtl1 <- eQTL_6mo
  #qtl2 <- eQTL_18mo

  # Make sure Ensembl Ids match up
  if(identical(qtl1$id, qtl2$id) == FALSE){
    stop("IDs do not match")
  }

  # combine and ask if peaks are local or distal
  QTLcompare <- qtl1 %>% mutate(
    AdditiveChr2 = qtl2$AdditiveChr,
    AdditivePos2 = qtl2$AdditivePos,
    AdditiveLOD2 = qtl2$AdditiveLOD,
    peak = AdditiveChr == AdditiveChr2
  )
  # Assign peak as "local"/"distant"
  QTLcompare[QTLcompare$peak == TRUE,]$peak <- "local"
  QTLcompare[QTLcompare$peak == FALSE,]$peak <- "distant"

  return(QTLcompare)
}

# Run comparison on all 3 combinations
# eQTL
eQTL_PeakComp6v18 <- CompQTL(eQTL_6mo, eQTL_18mo)
eQTL_PeakComp6v12 <- CompQTL(eQTL_6mo, eQTL_12mo)
eQTL_PeakComp12v18 <- CompQTL(eQTL_12mo, eQTL_18mo)
# pQTL
pQTL_PeakComp6v18 <- CompQTL(pQTL_6mo, pQTL_18mo)
pQTL_PeakComp6v12 <- CompQTL(pQTL_6mo, pQTL_12mo)
pQTL_PeakComp12v18 <- CompQTL(pQTL_12mo, pQTL_18mo)

# Create necessary objects for plotting map
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
map <- map_df_to_list(map = snps, pos_column = "bp")
K <- calc_kinship(probs, type = "loco", cores = 10)

# visualize
MakeEachQTL <- function(EnsID, symbol, type){

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
  plot_scan1(qtl_6mo, map)
  title(paste(symbol, "(", EnsID, ")", "@6mo"))
  p6 <- recordPlot()
  plot_scan1(qtl_12mo, map)
  title(paste(symbol, "(", EnsID, ")", "@12mo"))
  p12 <- recordPlot()
  plot_scan1(qtl_18mo, map)
  title(paste(symbol, "(", EnsID, ")", "@18mo"))
  p18 <- recordPlot()

  allplots <- list(p6,p12,p18)
  return(allplots)
  # see using:
  # allplots[[1]] for 6 months
  # allplots[[2]] for 12 months
  # allplots[[3]] for 18 months
}

# Pick outliers to look at individual eQTL per age maps---------------

# local 6v18
# high
pick1 <- eQTL_PeakComp6v18 %>%
  filter(peak == "local",
    AdditiveLOD2 > 35)
plot1 <- MakeEachQTL(pick1$id, pick1$symbol, "eQTL")
#show plots
plot1[[1]]
plot1[[2]]
plot1[[3]]
#low
pick2 <- eQTL_PeakComp6v18 %>%
  filter(peak == "local",
    AdditiveLOD2 < 10,
    AdditiveLOD > 23)
plot2 <- MakeEachQTL(pick2$id, pick2$symbol, "eQTL")
#show plots
plot2[[1]]
plot2[[2]]
plot2[[3]]

# distant 6v18
# high
pick1 <- eQTL_PeakComp6v18 %>%
  filter(peak == "distant",
    AdditiveLOD2 > 17)
plot1 <- MakeEachQTL(pick1$id, pick1$symbol, "eQTL")
#show plots
plot1[[1]]
plot1[[2]]
plot1[[3]]
#low
pick2 <- eQTL_PeakComp6v18 %>%
  filter(peak == "distant",
    AdditiveLOD > 21)
plot2 <- MakeEachQTL(pick2$id, pick2$symbol, "eQTL")
#show plots
plot2[[1]]
plot2[[2]]
plot2[[3]]

# local 6v12
# high
pick1 <- eQTL_PeakComp6v12 %>%
  filter(peak == "local",
    AdditiveLOD2 > 27,
    AdditiveLOD > 10,
    AdditiveLOD < 15)
plot1 <- MakeEachQTL(pick1$id, pick1$symbol, "eQTL")
#show plots
plot1[[1]]
plot1[[2]]
plot1[[3]]
#low
pick2 <- eQTL_PeakComp6v12 %>%
  filter(peak == "local",
    AdditiveLOD2 < 11,
    AdditiveLOD > 25)
plot2 <- MakeEachQTL(pick2$id, pick2$symbol, "eQTL")
#show plots
plot2[[1]]
plot2[[2]]
plot2[[3]]

# distant 6v12
# high
pick1 <- eQTL_PeakComp6v12 %>%
  filter(peak == "distant",
    AdditiveLOD2 > 18)
plot1 <- MakeEachQTL(pick1$id, pick1$symbol, "eQTL")
#show plots
plot1[[1]]
plot1[[2]]
plot1[[3]]
#low
pick2 <- eQTL_PeakComp6v12 %>%
  filter(peak == "distant",
    AdditiveLOD > 20)
plot2 <- MakeEachQTL(pick2$id, pick2$symbol, "eQTL")
#show plots
plot2[[1]]
plot2[[2]]
plot2[[3]]

# local 12v18
# high
pick1 <- eQTL_PeakComp12v18 %>%
  filter(peak == "local",
    AdditiveLOD2 > 35)
plot1 <- MakeEachQTL(pick1$id, pick1$symbol, "eQTL")
#show plots
plot1[[1]]
plot1[[2]]
plot1[[3]]
#low
pick2 <- eQTL_PeakComp12v18 %>%
  filter(peak == "local",
    AdditiveLOD2 < 15,
    AdditiveLOD > 28)
plot2 <- MakeEachQTL(pick2$id, pick2$symbol, "eQTL")
#show plots
plot2[[1]]
plot2[[2]]
plot2[[3]]

# distant 12v18
# high
pick1 <- eQTL_PeakComp12v18 %>%
  filter(peak == "distant",
    AdditiveLOD2 > 16.7)
plot1 <- MakeEachQTL(pick1$id, pick1$symbol, "eQTL")
#show plots
plot1[[1]]
plot1[[2]]
plot1[[3]]
#low
pick2 <- eQTL_PeakComp12v18 %>%
  filter(peak == "distant",
    AdditiveLOD > 17.5)
plot2 <- MakeEachQTL(pick2$id, pick2$symbol, "eQTL")
#show plots
plot2[[1]]
plot2[[2]]
plot2[[3]]

# pQTL--------------------------------------------------------------------------
# local 6v18
# high
pick1 <- pQTL_PeakComp6v18 %>%
  filter(peak == "local",
    AdditiveLOD2 > 31)
plot1 <- MakeEachQTL(pick1$id, pick1$symbol, "pQTL")
#show plots
plot1[[1]]
plot1[[2]]
plot1[[3]]
#low
pick2 <- pQTL_PeakComp6v18 %>%
  filter(peak == "local",
    AdditiveLOD > 34)
plot2 <- MakeEachQTL(pick2$id, pick2$symbol, "pQTL")
#show plots
plot2[[1]]
plot2[[2]]
plot2[[3]]

# distant 6v18
# high
pick1 <- pQTL_PeakComp6v18 %>%
  filter(peak == "distant",
    AdditiveLOD2 > 33.5)
plot1 <- MakeEachQTL(pick1$id, pick1$symbol, "pQTL")
#show plots
plot1[[1]]
plot1[[2]]
plot1[[3]]
#low
pick2 <- pQTL_PeakComp6v18 %>%
  filter(peak == "distant",
    AdditiveLOD > 32)
plot2 <- MakeEachQTL(pick2$id, pick2$symbol, "pQTL")
#show plots
plot2[[1]]
plot2[[2]]
plot2[[3]]

# local 6v12
# high
pick1 <- pQTL_PeakComp6v12 %>%
  filter(peak == "local",
    AdditiveLOD2 > 35)
plot1 <- MakeEachQTL(pick1$id, pick1$symbol, "pQTL")
#show plots
plot1[[1]]
plot1[[2]]
plot1[[3]]
#low
pick2 <- pQTL_PeakComp6v12 %>%
  filter(peak == "local",
    AdditiveLOD > 32,
    AdditiveLOD2 < 15)
plot2 <- MakeEachQTL(pick2$id, pick2$symbol, "pQTL")
#show plots
plot2[[1]]
plot2[[2]]
plot2[[3]]

# distant 6v12
# high
pick1 <- pQTL_PeakComp6v12 %>%
  filter(peak == "distant",
    AdditiveLOD2 > 49)
plot1 <- MakeEachQTL(pick1$id, pick1$symbol, "pQTL")
#show plots
plot1[[1]]
plot1[[2]]
plot1[[3]]
#low
pick2 <- pQTL_PeakComp6v12 %>%
  filter(peak == "distant",
    AdditiveLOD > 33)
plot2 <- MakeEachQTL(pick2$id, pick2$symbol, "pQTL")
#show plots
plot2[[1]]
plot2[[2]]
plot2[[3]]

# local 12v18
# high
pick1 <- pQTL_PeakComp12v18 %>%
  filter(peak == "local",
    AdditiveLOD2 > 27,
    AdditiveLOD < 15)
plot1 <- MakeEachQTL(pick1$id, pick1$symbol, "pQTL")
#show plots
plot1[[1]]
plot1[[2]]
plot1[[3]]
#low
pick2 <- pQTL_PeakComp12v18 %>%
  filter(peak == "local",
    AdditiveLOD > 35)
plot2 <- MakeEachQTL(pick2$id, pick2$symbol, "pQTL")
#show plots
plot2[[1]]
plot2[[2]]
plot2[[3]]

# distant 12v18
# high
pick1 <- pQTL_PeakComp12v18 %>%
  filter(peak == "distant",
    AdditiveLOD2 > 33.5)
plot1 <- MakeEachQTL(pick1$id, pick1$symbol, "pQTL")
#show plots
plot1[[1]]
plot1[[2]]
plot1[[3]]
#low
pick2 <- pQTL_PeakComp12v18 %>%
  filter(peak == "distant",
    AdditiveLOD > 48)
plot2 <- MakeEachQTL(pick2$id, pick2$symbol, "pQTL")
#show plots
plot2[[1]]
plot2[[2]]
plot2[[3]]
