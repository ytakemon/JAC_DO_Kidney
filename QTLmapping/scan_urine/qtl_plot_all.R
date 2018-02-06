#qsub -v script=qtl_plot_all Rsubmit_args.sh
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
  add_snp <- readRDS(paste0("./QTLscan/addscan_urine/Addsnps_", pheno, "_all.rds"))
  int_snp <- readRDS(paste0("./QTLscan/addscan_urine/Intsnps_", pheno, "_all.rds"))

  #188 set
  add_lod188 <- readRDS(paste0("./QTLscan/addscan_urine/Addscan_", pheno, "_188b.rds"))
  add_perm188 <- readRDS(paste0("./QTLscan/addscan_urine/Addperm_", pheno, "_188b.rds"))
  int_lod188 <- readRDS(paste0("./QTLscan/addscan_urine/Intscan_", pheno, "_188b.rds"))
  int_perm188 <- readRDS(paste0("./QTLscan/addscan_urine/Intperm_", pheno, "_188b.rds"))
  add_coef188 <- readRDS(paste0("./QTLscan/addscan_urine/Addcoef_", pheno, "_188b.rds"))
  int_coef188 <- readRDS(paste0("./QTLscan/addscan_urine/Intcoef_", pheno, "_188b.rds"))
  add_snp188 <- readRDS(paste0("./QTLscan/addscan_urine/Addsnps_", pheno, "_188b.rds"))
  int_snp188 <- readRDS(paste0("./QTLscan/addscan_urine/Intsnps_", pheno, "_188b.rds"))

  # Calc Age interactive factor
  Age_lod <- int_lod - add_lod
  Age_lod188 <- int_lod188 - add_lod188

  Age_perm <- add_perm - int_perm
  Age_perm188 <- add_perm188 - int_perm188

  #Make maps
  map_all <- map_df_to_list(map = MM_snps, pos_column = "pos")
  map_188 <- map_df_to_list(map = snps, pos_column = "pos")

  # QTL plot
  list <- c("add","int","add188", "int188", "Age", "Age188")
  # All add
  for (i in list){
    if (i == "add"){
      file <- paste0("./QTLscan/output/plots/Urine_", pheno, "_QTLmap_add.pdf")
      lod <- get("add_lod")
      perm <- get("add_perm")
      map <- get("map_all")
    } else if (i == "int"){
      file <- paste0("./QTLscan/output/plots/Urine_", pheno, "_QTLmap_int.pdf")
      lod <- get("int_lod")
      perm <- get("int_perm")
      map <- get("map_all")
    } else if (i == "add188"){
      file <- paste0("./QTLscan/output/plots/Urine_", pheno, "_QTLmap_add188b.pdf")
      lod <- get("add_lod188")
      perm <- get("add_perm188")
      map <- get("map_188")
    } else if (i == "int188"){
      file <- paste0("./QTLscan/output/plots/Urine_", pheno, "_QTLmap_int188b.pdf")
      lod <- get("int_lod188")
      perm <- get("int_perm188")
      map <- get("map_188")
    } else if (i == "add188"){
      file <- paste0("./QTLscan/output/plots/Urine_", pheno, "_QTLmap_add188b.pdf")
      lod <- get("add_lod188")
      perm <- get("add_perm188")
      map <- get("map_188")
    } else if (i == "Age"){
      file <- paste0("./QTLscan/output/plots/Urine_", pheno, "_QTLmap_Age.pdf")
      lod <- get("Age_lod")
      perm <- get("Age_perm")
      map <- get("map_all")
    } else if (i == "Age188"){
      file <- paste0("./QTLscan/output/plots/Urine_", pheno, "_QTLmap_Age188.pdf")
      lod <- get("Age_lod188")
      perm <- get("Age_perm188")
      map <- get("map_188")
    }
    print(paste("plotting QTL:", i))
    pdf(file = file, width = 12, height = 6)
    plot(lod, map)
    title(main = paste0(name, " QTL map (n = ", attributes(lod)$sample_size, ")"),
          sub = paste0("LOD threshold = ", signif(quantile(perm, 0.95)[1], digits = 3), " (0.05, 1000 permutations)"))
    abline( h = quantile(perm, 0.95)[1], col = "orange")
    dev.off()
  }

  # Coef plot
  # Get chr of highest peak
  add_chr <- max(add_lod, map_all)$chr
  int_chr <- max(int_lod, map_all)$chr
  add_chr188 <- max(add_lod188, map_188)$chr
  int_chr188 <- max(int_lod188, map_188)$chr
  Age_chr <- max(Age_lod, map_all)$chr
  Age_chr <- max(Age_lod, map_all)$chr

  list <- c("add","int","add188", "int188")
  for (i in list){
    if(i == "add"){
      lod <- get("add_lod")
      chr <- get("add_chr")
      map <- get("map_all")
      coef <- get("add_coef")
      file <- paste0("./QTLscan/output/plots/Urine_", pheno,  "_addQTLcoef_chr", chr,".pdf")
    } else if(i == "int"){
      lod <- get("int_lod")
      chr <- get("int_chr")
      map <- get("map_all")
      coef <- get("int_coef")
      file <- paste0("./QTLscan/output/plots/Urine_", pheno,  "_intQTLcoef_chr", chr,".pdf")
    } else if(i == "add188"){
      lod <- get("add_lod188")
      chr <- get("add_chr188")
      map <- get("map_188")
      coef <- get("add_coef188")
      file <- paste0("./QTLscan/output/plots/Urine_", pheno,  "_addQTLcoef188_chr", chr,".pdf")
    } else if(i == "int188"){
      lod <- get("int_lod188")
      chr <- get("int_chr188")
      map <- get("map_188")
      coef <- get("int_coef188")
      file <- paste0("./QTLscan/output/plots/Urine_", pheno,  "_intQTLcoef188_chr", chr,".pdf")
    }

    print(paste("Plotting allele coef: ", i))
    pdf(file = file, width = 12, height = 6)
    plot_coefCC(x = coef,
                map = map[chr],
                scan1_output=lod,
                legend = "topright")
    title(main = paste0(name, " QTL: Allele coefficient chr ", chr, " (n = ", attributes(lod)$sample_size, ")"))
    dev.off()
  }

  # Plot snps in peak lod interval
  # Link sql datebase
  query_genes <- create_gene_query_func("./qtl2_sqlite/mouse_genes_mgi.sqlite")

  list <- c("add","int","add188", "int188")
  for(i in list){
    if(i == "add"){
      lod <- get("add_lod")
      map <- get("map_all")
      getsnp <- get("add_snp")
      peak_chr <- max(lod, map)$chr
      peak_pos <- max(lod, map)$pos
      file <- paste0("./QTLscan/output/plots/Urine_", pheno,  "_addQTLgenes_chr", peak_chr,".pdf")
    } else if(i == "int"){
      lod <- get("int_lod")
      map <- get("map_all")
      getsnp <- get("int_snp")
      peak_chr <- max(lod, map)$chr
      peak_pos <- max(lod, map)$pos
      file <- paste0("./QTLscan/output/plots/Urine_", pheno,  "_intQTLgenes_chr", peak_chr,".pdf")
    } else if(i == "add188"){
      lod <- get("add_lod188")
      map <- get("map_188")
      getsnp <- get("add_snp188")
      peak_chr <- max(lod, map)$chr
      peak_pos <- (max(lod, map)$pos)/ 1e6 # need to convert to mega base
      file <- paste0("./QTLscan/output/plots/Urine_", pheno,  "_addQTLgenes188_chr", peak_chr,".pdf")
    } else if(i == "int188"){
      lod <- get("int_lod188")
      map <- get("map_188")
      getsnp <- get("int_snp188")
      peak_chr <- max(lod, map)$chr
      peak_pos <- (max(lod, map)$pos)/ 1e6  # need to convert to mega base
      file <- paste0("./QTLscan/output/plots/Urine_", pheno,  "_intQTLgenes188_chr", peak_chr,".pdf")
    }

    print(paste("Plotting genes under interval:", i))
    genes <- query_genes(peak_chr, peak_pos - 1, peak_pos + 1) # Pos has to be in megabases.

    pdf(file = file, width = 10, height = 6)
    par(mar=c(4.1, 4.1, 0.6, 0.6))
    plot(getsnp$lod, getsnp$snpinfo, drop_hilit=1.5, genes = genes)
    dev.off()
  }
}
