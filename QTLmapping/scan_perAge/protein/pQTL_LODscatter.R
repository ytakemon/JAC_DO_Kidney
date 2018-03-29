library(tidyverse)
library(grid)
library(ggsci)
options(dplyr.width = Inf) #override column limit
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
pQTL_best_6 <- read_csv("./QTLscan/output/6mo_pQTLBestperGene.csv")
pQTL_best_12 <- read_csv("./QTLscan/output/12mo_pQTLBestperGene.csv")
pQTL_best_18 <- read_csv("./QTLscan/output/18mo_pQTLBestperGene.csv")

# create function to generate df that made pQTL plots (incase)
# no threshold is set at the moment
Getgmb <- function(pQTL_best){
  chr_full <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X","Y", "MT")
  AddQTL <- pQTL_best %>%
              mutate(
                chr = factor(chr, levels = chr_full),
                AdditiveChr = factor(AdditiveChr, levels= chr_full)
              ) %>%
              filter(!(chr %in% c("Y","MT")))

  # Annotate postion with genes and save file for sharing
  save <- arrange(AddQTL, AdditiveChr, AdditivePos)
  # save annotated list for sharing
  #write.csv(save, "./QTLscan/output/Threshold6_pQTL_6mo.csv", row.names = FALSE, quote = FALSE)


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
pQTL_6mo <- Getgmb(pQTL_best_6)
pQTL_12mo <- Getgmb(pQTL_best_12)
pQTL_18mo <- Getgmb(pQTL_best_18)

# Compare two time points
CompQTL <- function(qtl1, qtl2){

  # remove from function once done
  #qtl1 <- pQTL_6mo
  #qtl2 <- pQTL_18mo

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
PeakComp6v18 <- CompQTL(pQTL_6mo, pQTL_18mo)
PeakComp6v12 <- CompQTL(pQTL_6mo, pQTL_12mo)
PeakComp12v18 <- CompQTL(pQTL_12mo, pQTL_18mo)

# Plot by LOD scores
# ? set threshold?

PlotbyLOD <- function(PeakComp, age1, age2, qtlType, peaktype){
  # remove after function is stable
  #PeakComp <- PeakComp6v18
  #peaktype options c("both","distant", "local")
  if(peaktype == "both"){
    plot <- ggplot(PeakComp, aes(x = AdditiveLOD, y = AdditiveLOD2, shape = peak, colour = peak))+
      geom_point(alpha = 0.3) +
      scale_colour_aaas() +
      geom_abline(intercept = 0, slope = 1, color="black") +
      labs( title = paste("LOD score comparision between", age1, "and", age2, "month", qtlType),
        subtitle = paste(nrow(PeakComp), "genes"))+
      scale_x_continuous(name = paste(age1,"month",qtlType,"LOD scores"))+
      scale_y_continuous(name = paste(age2,"month",qtlType,"LOD scores"))

  } else if(peaktype == "distant" | peaktype == "local"){
    plot <- PeakComp %>% filter(peak == peaktype) %>%
      ggplot(., aes(x = AdditiveLOD, y = AdditiveLOD2))+
      geom_point(alpha = 0.3) +
      scale_colour_aaas() +
      geom_abline(intercept = 0, slope = 1, color="black") +
      labs( title = paste(age1, "and", age2, "month", qtlType, "(",peaktype,"only)"),
        subtitle = paste(nrow(PeakComp %>% filter(peak == peaktype)), "genes"))+
      scale_x_continuous(name = paste(age1,"month",qtlType,"LOD scores"))+
      scale_y_continuous(name = paste(age2,"month",qtlType,"LOD scores"))
  } else {
    stop(paste("Do not recognist peak type", peaktype))
  }
  return(plot)
}

both6v18 <- PlotbyLOD(PeakComp = PeakComp6v18,
    age1 = "6",
    age2 = "18",
    qtlType = "pQTL",
    peaktype = "both")
both6v12 <- PlotbyLOD(PeakComp = PeakComp6v12,
    age1 = "6",
    age2 = "12",
    qtlType = "pQTL",
    peaktype = "both")
both12v18 <- PlotbyLOD(PeakComp = PeakComp12v18,
    age1 = "12",
    age2 = "18",
    qtlType = "pQTL",
    peaktype = "both")

local6v18 <- PlotbyLOD(PeakComp = PeakComp6v18,
    age1 = "6",
    age2 = "18",
    qtlType = "pQTL",
    peaktype = "local")

local6v12 <- PlotbyLOD(PeakComp = PeakComp6v12,
    age1 = "6",
    age2 = "12",
    qtlType = "pQTL",
    peaktype = "local")

local12v18 <- PlotbyLOD(PeakComp = PeakComp12v18,
    age1 = "12",
    age2 = "18",
    qtlType = "pQTL",
    peaktype = "local")

distant6v18 <- PlotbyLOD(PeakComp = PeakComp6v18,
    age1 = "6",
    age2 = "18",
    qtlType = "pQTL",
    peaktype = "distant")

distant6v12 <- PlotbyLOD(PeakComp = PeakComp6v12,
    age1 = "6",
    age2 = "12",
    qtlType = "pQTL",
    peaktype = "distant")

distant12v18 <- PlotbyLOD(PeakComp = PeakComp12v18,
    age1 = "12",
    age2 = "18",
    qtlType = "pQTL",
    peaktype = "distant")

pdf("./QTLscan/output/plots/pQTL_LODscore_comarison.pdf", width = 18, height = 12)
pushViewport(viewport(layout = grid.layout(12,18)))
print(both6v18, vp = viewport(layout.pos.row = 1:4, layout.pos.col = 1:8))
print(local6v18, vp = viewport(layout.pos.row = 1:4, layout.pos.col = 9:13))
print(distant6v18, vp = viewport(layout.pos.row = 1:4, layout.pos.col = 14:18))
print(both6v12, vp = viewport(layout.pos.row = 5:8, layout.pos.col = 1:8))
print(local6v12, vp = viewport(layout.pos.row = 5:8, layout.pos.col = 9:13))
print(distant6v12, vp = viewport(layout.pos.row = 5:8, layout.pos.col = 14:18))
print(both12v18, vp = viewport(layout.pos.row = 9:12, layout.pos.col = 1:8))
print(local12v18, vp = viewport(layout.pos.row = 9:12, layout.pos.col = 9:13))
print(distant12v18, vp = viewport(layout.pos.row = 9:12, layout.pos.col = 14:18))
dev.off()

CreateTable <- function(PeakComp, threshold, type){
  table <- PeakComp %>% filter(peak == type) %>%
    mutate(LODdiff = AdditiveLOD2 - AdditiveLOD) %>%
    filter(LODdiff > threshold | LODdiff < -threshold) %>%
    select(-t_gbm, -q_gbm)
  return(table)
}

# Create table
#local
Table6v18_local <- CreateTable(PeakComp6v18, 6, "local") %>%
  rename(AdditiveChr6mo = AdditiveChr,
    AdditiveLOD6mo = AdditiveLOD,
    AdditivePos6mo = AdditivePos,
    AdditiveChr18mo = AdditiveChr2,
    AdditiveLOD18mo = AdditiveLOD2,
    AdditivePos18mo = AdditivePos2)
Table6v12_local <- CreateTable(PeakComp6v12, 6, "local") %>%
  rename(AdditiveChr6mo = AdditiveChr,
    AdditiveLOD6mo = AdditiveLOD,
    AdditivePos6mo = AdditivePos,
    AdditiveChr12mo = AdditiveChr2,
    AdditiveLOD12mo = AdditiveLOD2,
    AdditivePos12mo = AdditivePos2)
Table12v18_local <- CreateTable(PeakComp12v18, 6, "local") %>%
  rename(AdditiveChr12mo = AdditiveChr,
    AdditiveLOD12mo = AdditiveLOD,
    AdditivePos12mo = AdditivePos,
    AdditiveChr18mo = AdditiveChr2,
    AdditiveLOD18mo = AdditiveLOD2,
    AdditivePos18mo = AdditivePos2)
#distant
Table6v18_distant <- CreateTable(PeakComp6v18, 6, "distant") %>%
  rename(AdditiveChr6mo = AdditiveChr,
    AdditiveLOD6mo = AdditiveLOD,
    AdditivePos6mo = AdditivePos,
    AdditiveChr18mo = AdditiveChr2,
    AdditiveLOD18mo = AdditiveLOD2,
    AdditivePos18mo = AdditivePos2)
Table6v12_distant <- CreateTable(PeakComp6v12, 6, "distant") %>%
  rename(AdditiveChr6mo = AdditiveChr,
    AdditiveLOD6mo = AdditiveLOD,
    AdditivePos6mo = AdditivePos,
    AdditiveChr12mo = AdditiveChr2,
    AdditiveLOD12mo = AdditiveLOD2,
    AdditivePos12mo = AdditivePos2)
Table12v18_distant <- CreateTable(PeakComp12v18, 6, "distant") %>%
  rename(AdditiveChr12mo = AdditiveChr,
    AdditiveLOD12mo = AdditiveLOD,
    AdditivePos12mo = AdditivePos,
    AdditiveChr18mo = AdditiveChr2,
    AdditiveLOD18mo = AdditiveLOD2,
    AdditivePos18mo = AdditivePos2)

# Write table
#local
write.csv(Table6v18_local, file = "./QTLscan/output/pQTL_LODdiff6_6v18_local.csv", row.names = FALSE)
write.csv(Table6v12_local, file = "./QTLscan/output/pQTL_LODdiff6_6v12_local.csv", row.names = FALSE)
write.csv(Table12v18_local, file = "./QTLscan/output/pQTL_LODdiff6_12v18_local.csv", row.names = FALSE)
#distant
write.csv(Table6v18_distant, file = "./QTLscan/output/pQTL_LODdiff6_6v18_distant.csv", row.names = FALSE)
write.csv(Table6v12_distant, file = "./QTLscan/output/pQTL_LODdiff6_6v12_distant.csv", row.names = FALSE)
write.csv(Table12v18_distant, file = "./QTLscan/output/pQTL_LODdiff6_12v18_distant.csv", row.names = FALSE)
