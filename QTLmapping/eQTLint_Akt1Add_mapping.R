setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
source("./Scripts/QTLmapping/miniDOQTL.R")
load("./RNAseq_data/DO188b_kidney.RData")
load("./shiny_annotation.RData")


# Get list of genes with trans eQTL
list <- read.csv("./QTLscan/output/Threshold8_eQTL_intAge.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == 12, ]
list <- list$symbol
level <- "mRNA"

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

for (i in 1:length(list)){
  gene <- annot.mrna[annot.mrna$symbol == list[i],]
  file_name <- paste0("./QTLscan/intscan_mrna_addAkt1/", gene$id, "_", gene$symbol, ".rds")
  file_name2 <- paste0("./QTLscan/intscan_mrna/Age/", gene$id, "_", gene$symbol, ".rds")
  file_name3 <- paste0("/hpcdata/cgd/QTL_mapping/kidney_combined/scanone_mrna/","ENSMUSG00000000001_Gnai3.rds")
  fit <- readRDS(paste0(file_name))
  fit2 <- readRDS(paste0(file_name2))
  fit3 <- readRDS(paste0(file_name3))

  # sign. thresholds based on 10000 permutations (100 permutations x 100 randomly selected gnes)
  sig.thr <- structure(c(6.28, 8.13, 9.14,
                       9.08, 11.32, 12.57), .Dim = c(3L,
                       2L), .Dimnames = list(c("0.63", "0.05", "0.01"), c("A", "X")))

  # plot
  plot.doqtl(fit, sig.thr = sig.thr, stat.name = "pheno1", main=gene, sig.col = c("yellow", "blue", "red"))
  plot.doqtl(fit2, sig.thr = sig.thr, stat.name = "pheno1", main=gene, sig.col = c("yellow", "blue", "red"))
}

plot.doqtl = function(x, stat.name = NULL,  sig.thr = NULL,
                      sig.col = "red", main = NULL, ...) {

  lod <- as.data.frame(x)
  lod$chr <- as.factor(substring(rownames(lod), 1, regexpr("_", rownames(lod)) - 1))
  lod$chr <- factor(lod$chr, levels(lod$chr)[c(1,12:19,2:11,20)])
  lod$mb <- as.numeric(substring(rownames(lod), regexpr("_", rownames(lod)) + 1))

  # Convert transcript and qtl position relative to chromosome positions
  # Convert to megabases
  chrlen <- sapply(split(lod$mb, lod$chr), max)
  chrsum <- c(0, cumsum(chrlen))
  names(chrsum) = names(chrlen)
  stat.name = stat.name

  # Get the call and arguments.
  call = match.call()

  # Get chr lengths and locations.  Create an x-axis based on Genome Mb.
  if(max(lod$mb, na.rm = TRUE) > 200) {
    lod$mb = lod$mb * 1e-6
  } # if(max(lod[,1]) > 200)

  mb = lod$mb
  gmb = mb
  unique.chr = as.character(unique(lod$chr))
  old.warn = options("warn")$warn
  options(warn = -1)
  unique.chr = unique.chr[order(as.numeric(unique.chr))]
  options(warn = old.warn)
  chrlen = chrlen[names(chrlen) %in% unique.chr]

  if(max(chrlen) > 200) {
    chrlen = chrlen * 1e-6
  } # if(max(chrlen) > 200)

  # Get the cumulative sum of the Chr lengths.
  chrsum = cumsum(chrlen)
  chrsum = c(0, chrsum)
  chrmid = chrsum[-length(chrsum)] + (diff(chrsum) * 0.5)

  # Add the preceeding chromosome lengths to each SNP position.
  for(c in 2:length(unique.chr)) {
    rows = which(lod$chr == unique.chr[c])
    gmb[rows] = gmb[rows] + chrsum[c]
  } # for(c)

  plot.column = which(colnames(lod) == stat.name)
  if(length(plot.column) == 0) {
    stop(paste("The stat.name of", stat.name, "was not found in the column names",
         "of the DOQTL object. Please verify that stat.name contains one of the",
         "column names in the DOQTL object."))
  } # if(length(plot.column) == 0)

  par(font = 2, font.lab = 2, font.axis = 2, las = 1, xaxs = "i",
      plt = c(0.12, 0.95, 0.05, 0.88))

  title = paste0(main$symbol, " Add: Akt1 (chr ", main$chr, ": ", main$start, "-", main$end,")")
  plot(gmb, lod[,1], col = 0, xlab = "", ylab = "lod", main = title,
       xaxt = "n", ylim = c(0, max(lod[,plot.column], na.rm = TRUE) * 1.05))

  usr = par("usr")
  rect(chrsum[2 * 1:(length(chrsum) * 0.5) - 1], usr[3],
       chrsum[2 * 1:(length(chrsum) * 0.5)], usr[4], col = rgb(0.8,0.8,0.8),
       border = NA)
  rect(usr[1], usr[3], usr[2], usr[4], border = 1)
  text(chrmid, 0.95 * usr[4], names(chrsum)[-1])

  lod = cbind(lod, gmb)
  lod = split(lod, lod$chr)
  lapply(lod, function(z) { points(z$gmb, z[,plot.column], type = "l", lwd = 2)})

  if(!is.null(sig.thr)) {
    add.sig.thr(sig.thr = sig.thr, sig.col = sig.col, chrsum = chrsum)
  }
} # plot.doqtl()
