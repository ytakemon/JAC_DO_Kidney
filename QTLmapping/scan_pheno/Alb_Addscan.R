library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO1045_kidney.Rdata")

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, MM_snps, pos_column = "marker")
MM_snps$chr <- as.character(MM_snps$chr)
MM_snps$chr[MM_snps$chr=="X"] <- "20"
map <- map_df_to_list(map = MM_snps, pos_column = "marker")

addcovar <- model.matrix(~ Sex + Generation + Cohort + DOB, data = samples)

# All phenotype QTL scans will be put off until current paper is published




lod <- scan1(genoprobs=probs,
             kinship=Glist2,
             pheno=expr.protein[present,p],
             addcovar=addcovar[,-1],
             cores=10, reml=TRUE)


for (p in plist) {

  cat("Scanning ",which(p==plist)," out of ",length(plist),"\n")

  present <- !is.na(expr.protein[,p])
  for (j in 1:20)
    Glist2 <- Glist[[j]][present,present]
  if (length(unique(annot.samples$Generation[present]))>1){
    addcovar <- model.matrix(~ Sex + Age + Generation + Protein.Batch + Protein.Channel, data=annot.samples[present, ])
  } else {
    addcovar <- model.matrix(~ Sex + Age + Protein.Batch + Protein.Channel, data=annot.samples[present, ])
  }

  # lod score
  lod <- scan1(genoprobs=probs[present,],
               kinship=Glist2,
               pheno=expr.protein[present,p],
               addcovar=addcovar[,-1],
               cores=10, reml=TRUE)

  # save lod object
  file_name <- paste0("./QTLscan/addscan_prot_pbatch/", annot.protein$id[p], "_", annot.protein$symbol[p], ".rds")
  saveRDS(lod, file=file_name)
}


operm <- scan1perm(pr, iron$pheno, Xcovar=Xcovar, n_perm=1000)
operm <- scan1perm(pr, iron$pheno, Xcovar=Xcovar, n_perm=1000, cores=0)

probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="X"] <- "20"
map <- map_df_to_list(map = snps, pos_column = "bp")


plist <- plist[plist<=ncol(expr.protein)]




create.plot <- function(gene.name, level, plotType, chr) {
  # if not ready, plot nothing
  if (is.null(gene.name) | is.null(level) | is.null(plotType) | is.null(chr)) return(NULL)

  file.location <- file.path(gene.name, level, plotType)
  if (!file.exists(file.location)) {nodata(); return(NULL)}
  fit <- readRDS(file.location)

  title <- gene.name
  if (level=="protein") title <- toupper(title)

  # lod plots
  if (plotType == 0 | plotType == 1 | plotType == 2) {

    # sign. thresholds based on 10000 permutations (100 permutations x 100 randomly selected gnes)
    if (level=="protein") {
      sig.thr <- structure(c(6.33, 8.16, 9.22,
                           9.22, 11.45, 12.58), .Dim = c(3L,
                           2L), .Dimnames = list(c("0.63", "0.05", "0.01"), c("A", "X")))
    } else {
      sig.thr <- structure(c(6.28, 8.13, 9.14,
                           9.08, 11.32, 12.57), .Dim = c(3L,
                           2L), .Dimnames = list(c("0.63", "0.05", "0.01"), c("A", "X")))
    }

    if (chr == "all") {
      plot.doqtl(fit, sig.thr = sig.thr, main=title, sig.col = c("yellow", "blue", "red"))
    } else {
      if (chr!="X") {
        coefplot(fit, chr, main=title)
      } else {
        coefplot(fit, chr, main=paste(title, "females"), sex="F")
      }
    }
  }
