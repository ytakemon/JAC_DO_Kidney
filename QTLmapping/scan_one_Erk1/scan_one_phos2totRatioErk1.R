# pbsnodes -a
# qsub -q short -X -l nodes=cadillac079:ppn=3,walltime=3:59:00 -I
# qsub -v script=scan_one_phos2totRatioErk1 Rsubmit_args.sh

library(qtl2)
library(qtl2convert)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO1045_kidney.Rdata")
erk1 <- read.delim("./Phenotype/phenotypes/JAC_WB_kidney_ERK.txt")
pheno <- "Phospho_ERK1_ratio"

# Cleanup data and subset to match samples
# match pheno to samples
erk1 <- erk1[complete.cases(erk1[,pheno]) & complete.cases(erk1[,"Phospho_ERK1"]),]
erk1 <- erk1[erk1$ID %in% samples$Mouse.ID, ]
erk1$duplicated <- (duplicated(erk1$ID) | duplicated(erk1$ID, fromLast = TRUE))
erk1 <- erk1[erk1$duplicated == FALSE,]
rownames(erk1) <- erk1$ID

# match samples to pheno
sub_samples <- samples[rownames(samples) %in% rownames(erk1),]

#match genoprobs to samples
sub_genoprobs <- genoprobs[rownames(genoprobs) %in% rownames(sub_samples),,]

# Check identical
identical(rownames(erk1), rownames(sub_samples))
identical(rownames(erk1), rownames(sub_genoprobs))

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(sub_genoprobs, MM_snps, pos_column = "pos")
K <- calc_kinship(probs, "loco", cores=3)
MM_snps$chr <- as.character(MM_snps$chr)
MM_snps$chr[MM_snps$chr=="20"] <- "X"
snps <- MM_snps[dimnames(sub_genoprobs)[[3]],]
map <- map_df_to_list(map = snps, pos_column = "pos")

addcovar <- model.matrix(~ Sex + Generation + Cohort.Age.mo, data = sub_samples)

lod <- scan1(genoprobs=probs,
             kinship=K,
             pheno=log(erk1[,pheno, drop = FALSE]),
             addcovar=addcovar[,-1],
             cores=3, reml=TRUE)

# save lod
file_name <- paste0("./QTLscan/addscan_phenotype/", pheno, ".rds")
saveRDS(lod, file=file_name)

# plot
# load lod and perms
pheno <- "Phospho_ERK1_ratio"
file_name <- paste0("./QTLscan/addscan_phenotype/", pheno, ".rds")
lod <- readRDS(file_name)
file_name <- paste0("./QTLscan/addscan_phenotype/", pheno, "_perm.rds")
perm <- readRDS(file_name)

pdf(paste0("./QTLscan/output/plots/", pheno,  "_qtl_map.pdf"), width = 12, height = 6)
plot(lod, map)
title(main = paste0("Phospho-Erk1 /  QTL map"),
      sub = paste0("LOD threshold = ", signif(summary(perm)[1], digits = 3), " (0.05, 1000 permutations)"))
abline( h = summary(perm)[1], col = "orange")
dev.off()
