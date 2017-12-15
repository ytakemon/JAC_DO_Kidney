# pbsnodes -a
# qsub -q short -X -l nodes=cadillac031:ppn=3,walltime=3:59:00 -I
# qsub -v script=scan_one_Mapk3 Rsubmit_args.sh

library(qtl2)
library(qtl2convert)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")
pheno <- "MAPK3_proteome"

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "pos")
K <- calc_kinship(probs, "loco", cores=3)
snps$chr <- as.character(snps$chr)
snps$chr[snps$chr=="20"] <- "X"
map <- map_df_to_list(map = snps, pos_column = "pos")

addcovar <- model.matrix(~ Sex + Generation + Age + Protein.Batch + Protein.Channel , data = annot.samples)

prot <- "Mapk3"
prot <- annot.protein[annot.protein$symbol == prot,]

lod <- scan1(genoprobs=probs,
             kinship=K,
             pheno=expr.protein[,prot$id],
             addcovar=addcovar[,-1],
             cores=3, reml=TRUE)

# save lod
file_name <- paste0("./QTLscan/addscan_phenotype/", pheno, ".rds")
saveRDS(lod, file=file_name)

# plot
# load lod and perms
pheno <- "MAPK3_proteome"
file_name <- paste0("./QTLscan/addscan_phenotype/", pheno, ".rds")
lod <- readRDS(file_name)
file_name <- paste0("./QTLscan/addscan_phenotype/", pheno, "_perm.rds")
perm <- readRDS(file_name)

pdf(paste0("./QTLscan/output/plots/", pheno, "_qtl_map.pdf"), width = 12, height = 6)
plot(lod, map)
title(main = paste0("MAPK3 (proteomics) QTL map"),
      sub = paste0("LOD threshold = ", signif(summary(perm)[1], digits = 3), " (0.05, 1000 permutations)"))
abline( h = summary(perm)[1], col = "orange")
dev.off()
