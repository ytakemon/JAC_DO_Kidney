library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(qtl2)
library(dplyr)
library(gridExtra)
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
  add_coef <- readRDS(paste0("./QTLscan/addscan_urine/Addcoef_", pheno, "_all.rds"))
  int_coef <- readRDS(paste0("./QTLscan/addscan_urine/Intcoef_", pheno, "_all.rds"))

  #188 set
  add_lod188 <- readRDS(paste0("./QTLscan/addscan_urine/Addscan_", pheno, "_188b.rds"))
  add_perm188 <- readRDS(paste0("./QTLscan/addscan_urine/Addperm_", pheno, "_188b.rds"))
  int_lod188 <- readRDS(paste0("./QTLscan/addscan_urine/Intscan_", pheno, "_188b.rds"))
  int_perm188 <- readRDS(paste0("./QTLscan/addscan_urine/Intperm_", pheno, "_188b.rds"))
  add_coef188 <- readRDS(paste0("./QTLscan/addscan_urine/Addcoef_", pheno, "_188b.rds"))
  int_coef188 <- readRDS(paste0("./QTLscan/addscan_urine/Intcoef_", pheno, "_188b.rds"))

  #Make maps
  map_all <- map_df_to_list(map = MM_snps, pos_column = "pos")
  map_188 <- map_df_to_list(map = snps, pos_column = "pos")

  # QTL plot
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

  # Coef plot
  # Get chr of highest peak
  add_chr <- max(add_lod, map_all)$chr
  int_chr <- max(int_lod, map_all)$chr
  add_chr188 <- max(add_lod188, map_188)$chr
  int_chr188 <- max(int_lod188, map_188)$chr

  #plot
  pdf(paste0("./QTLscan/output/plots/Urine_", pheno,  "_addQTLcoef_chr", add_chr,".pdf"), width = 12, height = 6)
  plot_coefCC(x = add_coef,
              map = map_all[add_chr],
              scan1_output=add_lod)
  legend("topright", col=CCcolors, names(CCcolors), ncol=2, lwd=2, bg="gray95")
  title(main = paste0(name, " QTL: Allele coefficient chr ", add_chr, " (n = ", attributes(add_lod)$sample_size, ")"))
  dev.off()

  pdf(paste0("./QTLscan/output/plots/Urine_", pheno,  "_intQTLcoef_chr", int_chr,".pdf"), width = 12, height = 6)
  plot_coefCC(x = int_coef,
              map = map_all[int_chr],
              scan1_output=int_lod)
  legend("topright", col=CCcolors, names(CCcolors), ncol=2, lwd=2, bg="gray95")
  title(main = paste0(name, " QTL: Allele coefficient chr ", int_chr, " (n = ", attributes(int_lod)$sample_size, ")"))
  dev.off()

  pdf(paste0("./QTLscan/output/plots/Urine_", pheno,  "_addQTLcoef188_chr", add_chr188,".pdf"), width = 12, height = 6)
  plot_coefCC(x = add_coef188,
              map = map_188[add_chr188],
              scan1_output=add_lod188)
  legend("topright", col=CCcolors, names(CCcolors), ncol=2, lwd=2, bg="gray95")
  title(main = paste0(name, " QTL: Allele coefficient chr ", add_chr188, " (n = ", attributes(add_lod188)$sample_size, ")"))
  dev.off()

  pdf(paste0("./QTLscan/output/plots/Urine_", pheno,  "_intQTLcoef188_chr", int_chr188,".pdf"), width = 12, height = 6)
  plot_coefCC(x = int_coef188,
              map = map_188[int_chr188],
              scan1_output=int_lod188)
  legend("topright", col=CCcolors, names(CCcolors), ncol=2, lwd=2, bg="gray95")
  title(main = paste0(name, " QTL: Allele coefficient chr ", int_chr188, " (n = ", attributes(int_lod188)$sample_size, ")"))
  dev.off()
}
