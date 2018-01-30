library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(qtl2)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO1045_kidney.Rdata")
load("./RNAseq_data/DO188b_kidney.RData")

pheno_list <- c("alb", "phs")

for (pheno in pheno_list){
  # Assign name
  if (pheno == "alb"){
    name <- "Albumin"
  } else {
    name <- "Phosphate"
  }

  # Load data:
  #All
  add_lod <- readRDS(paste0("./QTLscan/addscan_urine/Addscan_", pheno, "_all.rds"))
  add_perm <- readRDS(paste0("./QTLscan/addscan_urine/Addperm_", pheno, "_all.rds"))
  int_lod <- readRDS(paste0("./QTLscan/addscan_urine/Intscan_", pheno, "_all.rds"))
  int_perm <- readRDS(paste0("./QTLscan/addscan_urine/Intperm_", pheno, "_all.rds"))

  #188 set
  add_lod188 <- readRDS(paste0("./QTLscan/addscan_urine/Addscan_", pheno, "_188b.rds"))
  add_perm188 <- readRDS(paste0("./QTLscan/addscan_urine/Addperm_", pheno, "_188b.rds"))
  int_lod188 <- readRDS(paste0("./QTLscan/addscan_urine/Intscan_", pheno, "_188b.rds"))
  int_perm188 <- readRDS(paste0("./QTLscan/addscan_urine/Intperm_", pheno, "_188b.rds"))

  #Make maps
  map_all <- map_df_to_list(map = map, pos_column = "pos")
  map_188 <- map_df_to_list(map = snps, pos_column = "pos")

  #plot
  pdf(paste0("./QTLscan/output/plots/Urine_", pheno, "_QTLmap_add.pdf"), width = 12, height = 6)
  plot(add_lod, map_all)
  title(main = paste0(name, " QTL map (n = ", attributes(add_lod)$sample_size, ")"),
        sub = paste0("LOD threshold = ", signif(quantile(add_perm, 0.95)[1], digits = 3), " (0.05, 1000 permutations)"))
  abline( h = quantile(add_perm, 0.95)[1], col = "orange")
  dev.off()

  pdf(paste0("./QTLscan/output/plots/Urine_", pheno, "_QTLmap_int.pdf"), width = 12, height = 6)
  plot(int_lod, map_all)
  title(main = paste0(name, " QTL map (n = ", attributes(int_lod)$sample_size, ")"),
        sub = paste0("LOD threshold = ", signif(quantile(int_perm, 0.95)[1], digits = 3), " (0.05, 1000 permutations)"))
  abline( h = quantile(int_perm, 0.95)[1], col = "orange")
  dev.off()

  pdf(paste0("./QTLscan/output/plots/Urine_", pheno, "_QTLmap_add188b.pdf"), width = 12, height = 6)
  plot(add_lod188, map_188)
  title(main = paste0(name, " QTL map (n = ", attributes(add_lod188)$sample_size, ")"),
        sub = paste0("LOD threshold = ", signif(quantile(add_perm, 0.95)[1], digits = 3), " (0.05, 1000 permutations)"))
  abline( h = quantile(add_perm188, 0.95)[1], col = "orange")
  dev.off()

  pdf(paste0("./QTLscan/output/plots/Urine_", pheno, "_QTLmap_int188b.pdf"), width = 12, height = 6)
  plot(int_lod188, map_188)
  title(main = paste0(name, " QTL map (n = ", attributes(int_lod188)$sample_size, ")"),
        sub = paste0("LOD threshold = ", signif(quantile(int_perm188, 0.95)[1], digits = 3), " (0.05, 1000 permutations)"))
  abline( h = quantile(int_perm, 0.95)[1], col = "orange")
  dev.off()
}






quantile(add_perm, 0.90)
quantile(add_perm, 0.65)




# plot
# load lod188 and perms
pheno <- "Total_ERK1"
file_name <- paste0("./QTLscan/addscan_phenotype/", pheno, ".rds")
lod <- readRDS(file_name)
file_name <- paste0("./QTLscan/addscan_phenotype/", pheno, "_perm.rds")
perm <- readRDS(file_name)

pdf(paste0("./QTLscan/output/plots/", pheno, "_qtl_map.pdf"), width = 12, height = 6)
plot(lod, map)
title(main = paste0(name, " QTL map (n = ", attribute(add_lod), ")"),
      sub = paste0("LOD threshold = ", signif(summary(perm)[1], digits = 3), " (0.05, 1000 permutations)"))
abline( h = summary(perm)[1], col = "orange")
dev.off()
