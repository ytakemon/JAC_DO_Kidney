# Alb_
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO1045_kidney.Rdata")

# Subset pheontype: 6mo alb
pheno <- Upheno[Upheno$study == "Cross-sectional",]
pheno <- pheno[-1,]
pheno_cr <- pheno[,c("Mouse.ID", "cr.u.6", "cr.u.12", "cr.u.18")]
pheno_ma <- pheno[,c("Mouse.ID", "ma.u.6", "ma.u.12", "ma.u.18")]

# colapse to one column and combine
pheno_cr <- pheno_cr %>% mutate(
              cr.u.all = coalesce(cr.u.6, cr.u.12, cr.u.18)
            )
pheno_ma <- pheno_ma %>% mutate(
              ma.u.all = coalesce(ma.u.6, ma.u.12, ma.u.18)
            )

pheno <- pheno_cr[,c("Mouse.ID","cr.u.all")]
pheno$ma.u.all <- pheno_ma$ma.u.all

# filter out samples with missing phenotype
pheno <- pheno[!is.na(pheno[,2] & pheno[,3]),]






samples <- samples[samples$Cohort == "Cross-Sectional",]
samples$Mouse.ID <- as.character(samples$Mouse.ID)

pheno$Mouse.ID[!(pheno$Mouse.ID %in% samples$Mouse.ID)]


pheno6 <- Upheno[,c("Mouse.ID", "cr.u.6", "ma.u.6", "study")]
present <- !is.na(pheno6[,2] & pheno6[,3])
pheno6 <- pheno6[present,]


# Subset dataset
genoprobs <- genoprobs[pheno6$Mouse.ID,,]
samples <- samples[pheno6$Mouse.ID,]

# prepare data for qtl2
MM_snps$chr <- as.character(MM_snps$chr)
probs <- probs_doqtl_to_qtl2(genoprobs, MM_snps, pos_column = "pos")
K <- calc_kinship(probs, type = "loco")
map <- map_df_to_list(map = MM_snps, pos_column = "pos")
samples$cr <- pheno6[,2]
addcovar <- model.matrix(~ Sex + Generation + Cohort + DOB, data = samples)

# scan
lod <- scan1(genoprobs=probs,
             kinship=K,
             pheno=pheno6[,3],
             addcovar=addcovar[,-1],
             cores=20,
             reml=TRUE)


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
