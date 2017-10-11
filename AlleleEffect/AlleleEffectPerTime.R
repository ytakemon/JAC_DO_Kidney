library(shiny)
library(qtl)
library(qtlcharts)
library(ggplot2)
source("/hpcdata/cgd/shinyapps/kidney/miniDOQTL.R")
load("/hpcdata/cgd/shinyapps/kidney/shiny_annotation.RData")

plotType <- 0
gene.name <- "Slc6a20b"
level <- "mRNA"
chr <- 9

# for given MGI symbol, find Ensembl ids
other.ids <- function(gene.name, level) {
  if (level == "mRNA") {
    sel <- which(mRNA.list$symbol == gene.name)[1]
    if (!is.na(sel)) return(mRNA.list[sel,]) else return(c(NA,NA,NA))
  }
  if (level == "protein") {
    sel <- which(protein.list$symbol == gene.name)[1]
    if (!is.na(sel)) return(protein.list[sel,]) else return(c(NA,NA,NA))
  }
}

m.file.location <- paste0("/hpcdata/cgd/QTL_mapping/kidney_combined/scanone_mrna/", other.ids(gene.name, "mRNA")[[1]],"_",gene.name,".rds")
p.file.location <- paste0("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/QTLscan/addscan_mrna/", other.ids(gene.name, "mRNA")[[1]],"_",gene.name,".rds")

if (file.exists(m.file.location)){
  m.fit <- readRDS(m.file.location)
}
if (file.exists(p.file.location)){
  p.fit <- readRDS(p.file.location)
}

title <- gene.name
if (level=="protein") title <- toupper(title)

if (level=="protein") {
  sig.thr <- structure(c(6.33, 8.16, 9.22,
                       9.22, 11.45, 12.58), .Dim = c(3L,
                       2L), .Dimnames = list(c("0.63", "0.05", "0.01"), c("A", "X")))
} else {
  sig.thr <- structure(c(6.28, 8.13, 9.14,
                       9.08, 11.32, 12.57), .Dim = c(3L,
                       2L), .Dimnames = list(c("0.63", "0.05", "0.01"), c("A", "X")))
}

coefplot(m.fit, chr, main=paste(title, "females"), sex="F")
