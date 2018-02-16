#qsub -v script=qtl_plot_all Rsubmit_args.sh
library(qtl2convert)
library(qtl2)
library(dplyr)
library(gridExtra)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO1045_kidney.Rdata")

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

  if(pheno == "alb"){
    int_snp1 <- readRDS(paste0("./QTLscan/addscan_urine/Intsnps_", pheno, "_all2.rds"))
    int_snp2 <- readRDS(paste0("./QTLscan/addscan_urine/Intsnps_", pheno, "_all3.rds"))
  }

  # Calc Age interactive factor
  Age_lod <- int_lod - add_lod
  Age_perm <- sort(int_perm) - sort(add_perm)

  #Make maps
  map_all <- map_df_to_list(map = MM_snps, pos_column = "pos")

  # Link sql datebase
  query_genes <- create_gene_query_func("./qtl2_sqlite/mouse_genes_mgi.sqlite")

  # Get chr of highest peak
  add_chr <- max(add_lod, map_all)$chr
  int_chr <- max(int_lod, map_all)$chr
  Age_chr <- max(Age_lod, map_all)$chr

  # QTL plot
  list <- c("add", "int", "Age")
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
    } else if (i == "Age"){
      file <- paste0("./QTLscan/output/plots/Urine_", pheno, "_QTLmap_Age.pdf")
      lod <- get("Age_lod")
      perm <- get("Age_perm")
      map <- get("map_all")
    }
    print(paste("plotting QTL map:", i))
    pdf(file = file, width = 12, height = 6)
    plot(lod, map)
    title(main = paste0(name," ", i,"QTL map (n = ", attributes(lod)$sample_size, ")"),
          sub = paste0("LOD threshold = ", signif(quantile(perm, 0.95)[1], digits = 3), " (0.05, 1000 permutations)"))
    abline( h = quantile(perm, 0.95)[1], col = "orange")
    dev.off()
  }

  # Coef plot
  list <- c("add","int")
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
    }

    print(paste("Plotting allele coef: ", i))
    pdf(file = file, width = 12, height = 6)
    plot_coefCC(x = coef,
                map = map[chr],
                scan1_output=lod,
                legend = "topleft")
    title(main = paste0(name," ", i, "QTL: Allele coefficient chr ", chr, " (n = ", attributes(lod)$sample_size, ")"))
    dev.off()
  }

  # Plot snps in peak lod interval
  list <- c("add","int")
  for(i in list){
    if(i == "add"){
      lod <- get("add_lod")
      map <- get("map_all")
      getsnp <- get("add_snp")
      peak_chr <- get("add_chr")
      bayesint <- bayes_int(lod, map, peak_chr)
      start <- bayesint[,"ci_lo"] # start and end must match or be within initial scan1snp interval
      end <- bayesint[,"ci_hi"]
      file <- paste0("./QTLscan/output/plots/Urine_", pheno,  "_addQTLgenes_chr", peak_chr,".pdf")
    } else if(i == "int"){
      lod <- get("int_lod")
      map <- get("map_all")
      getsnp <- get("int_snp")
      peak_chr <- get("int_chr")
      bayesint <- bayes_int(lod, map, peak_chr)
      start <- bayesint[,"ci_lo"]

      # special case for int phs
      if (pheno == "phs"){
        end <- bayesint[,"pos"]+10
      } else {
        end <- bayesint[,"pos"]
      }
      file <- paste0("./QTLscan/output/plots/Urine_", pheno,  "_intQTLgenes_chr", peak_chr,".pdf")

      # special case for int albumin for closer look at gene under peak
      if (pheno == "alb" & i == "int"){
        getsnp1 <- get("int_snp1")
        start1 <- 15
        end1 <- 17
        file1 <- paste0("./QTLscan/output/plots/Urine_", pheno,  "_intQTLgenes_chr", peak_chr,"_mag1.pdf")

        getsnp2 <- get("int_snp2")
        start2 <- 24
        end2 <- 26
        file2 <- paste0("./QTLscan/output/plots/Urine_", pheno,  "_intQTLgenes_chr", peak_chr,"_mag2.pdf")
      }
    }

    print(paste("Plotting genes under interval:", i))
    genes <- query_genes(peak_chr, start, end) # Pos has to be in megabases.
    pdf(file = file, width = 10, height = 6)
    par(mar=c(4.1, 4.1, 1.5, 0.6))
    plot(getsnp$lod, getsnp$snpinfo, drop_hilit=1.5, genes = genes)
    title(main = paste0(name," ", i, "QTL: Genes under chr ", chr, " (Peak pos: ", bayesint[,"pos"],")"))
    dev.off()

    if (pheno == "alb" & i == "int"){
      print(paste("Plotting genes under interval:", i, "1"))
      genes <- query_genes(peak_chr, start1, end1) # Pos has to be in megabases.
      pdf(file = file1, width = 10, height = 6)
      par(mar=c(4.1, 4.1, 1.5, 0.6))
      plot(getsnp1$lod, getsnp1$snpinfo, drop_hilit=1.5, genes = genes)
      title(main = paste0(name," ", i, "QTL: Genes under chr ", chr, " Mag1 (Peak pos: ", bayesint[,"pos"],")"))
      dev.off()

      print(paste("Plotting genes under interval:", i, "2"))
      genes <- query_genes(peak_chr, start2, end2) # Pos has to be in megabases.
      pdf(file = file2, width = 10, height = 6)
      par(mar=c(4.1, 4.1, 1.5, 0.6))
      plot(getsnp2$lod, getsnp2$snpinfo, drop_hilit=1.5, genes = genes)
      title(main = paste0(name," ", i, "QTL: Genes under chr ", chr, " Mag2 (Peak pos: ", bayesint[,"pos"],")"))
      dev.off()
    }
  }
}
