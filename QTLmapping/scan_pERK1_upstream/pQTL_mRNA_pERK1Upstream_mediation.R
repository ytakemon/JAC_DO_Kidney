# qsub -v script=pQTL_mRNA_pERK1Upstream_mediation Rsubmit_args.sh

# R/3.4.1
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")
library(ggplot2)
library(dplyr)
library(scales)

# Get list of genes with trans pQTL
list <- read.csv("./QTLscan/output/Threshold6_pQTL_intAge_pbatch.csv", header = TRUE, stringsAsFactors = FALSE)
chr <- "7"
list <- list[list$IntAgeChr == chr, ]

# Get list of query genes
query_list <- read.csv("./SpecificQ/pERK1_addtotERK1_LODpeak_genes.csv", header = TRUE, stringsAsFactors = FALSE)
E_query_list <- query_list[query_list$ensembl_gene_id %in% annot.mrna$id,]
P_query_list <- query_list[query_list$ensembl_gene_id %in% annot.protein$gene_id,]

for (g in 1:nrow(E_query_list)){

  # Get candidate mediator
  med_query <- E_query_list$ensembl_gene_id[g]
  med_genename <- E_query_list$mgi_symbol[g]

  # parameters
  addscan.dir <- paste0("./QTLscan/addscan_mrna_pERK1Upstream/")
  intscan.dir.Age <-  paste0("./QTLscan/intscan_mrna_pERK1Upstream/")
  output.file <- paste0("./QTLscan/output/pQTLBestperGene_",med_query,"_","_mrna_pERK1Upstream_thr6_chr7.csv")

  annot.protein <- annot.protein[annot.protein$id %in% list$id,]
  output <- annot.protein[,c(1:6,10)]
  output$AdditiveLOD <- output$AdditivePos <-  output$AdditiveChr <- NA
  output$IntAgeLODFull <- output$IntAgeLODDiff <- output$IntAgePos <- output$IntAgeChr <- NA

  file.name <- function(i) paste0("Add_",med_query,"_Medscan_",output$id[i],"_",output$symbol[i],".rds")

  for (i in 1:nrow(output)) {

    if (i %% 10 == 0) print(i)

    if(file.exists(paste0(addscan.dir,file.name(i))) && file.exists(paste0(intscan.dir.Age, file.name(i)))){
      fit0 <- readRDS(paste0(addscan.dir,file.name(i)))
      fitAge <- readRDS(paste0(intscan.dir.Age, file.name(i)))
    }else{
      next
    }

    # additive scan
    dt <- data.frame(AdditiveChr=sapply(strsplit(rownames(fit0),"_"), "[", 1),
                     AdditivePos=as.numeric(sapply(strsplit(rownames(fit0),"_"), "[", 2)),
                     AdditiveLOD=fit0[,1], stringsAsFactors=FALSE)

    # find max lod score of chr 7 ONLY!
    dt2 <- dt[dt$AdditiveChr == chr,]
    dt2 <- dt %>% group_by(AdditiveChr) %>%
          summarize(AdditivePos = AdditivePos[which.max(AdditiveLOD)[1]],
                     AdditiveLOD = max(AdditiveLOD)) %>%
           arrange(-AdditiveLOD)
    output[i, c("AdditiveChr", "AdditivePos", "AdditiveLOD")] <- dt2[dt2$AdditiveChr == chr,]

    # int. scan - Age
    stopifnot(rownames(fit0) == rownames(fitAge))
    dt <- data.frame(IntAgeChr=sapply(strsplit(rownames(fit0),"_"), "[", 1),
                     IntAgePos=as.numeric(sapply(strsplit(rownames(fit0),"_"), "[", 2)),
                     IntAgeLODDiff=fitAge[,1] - fit0[,1],
                     IntAgeLODFull = fitAge[,1],
                     stringsAsFactors=FALSE)

    # find max lod score of chr7 ONLY!
    dt2 <- dt %>% group_by(IntAgeChr) %>%
           summarize(IntAgePos = IntAgePos[which.max(IntAgeLODFull)[1]],
                     IntAgeLODDiff=IntAgeLODDiff[which.max(IntAgeLODFull)[1]],
                     IntAgeLODFull = max(IntAgeLODFull)) %>%
           arrange(-IntAgeLODDiff)
    output[i, c("IntAgeChr", "IntAgePos", "IntAgeLODDiff", "IntAgeLODFull")] <- dt2[dt2$IntAgeChr == chr,]

      # collect rows into one data frame
      write.csv(output, file=output.file, row.names=FALSE)

      # compare & make plot:
      list <- arrange(list, id)
      list_add <- arrange(output, id)

      if (identical(list$id, list_add$id)){
        compare <- list[,colnames(list) %in% c("id", "gene_id", "symbol", "IntAgeChr", "IntAgePos", "IntAgeLODDiff")]
        compare$addIntAgeChr <- list_add$IntAgeChr
        compare$addIntAgePos <- list_add$IntAgePos
        compare$addIntAgeLODDiff <- list_add$IntAgeLODDiff
        compare <- compare[complete.cases(compare$addIntAgeChr),]
        write.csv(compare, file=paste0("./QTLscan/output/pERK1_upstream/mRNA/pERK1_upstream_mRNA_mediator_", med_query ,"_", med_genename, "_compare_chr7_thr6.csv"), row.names = FALSE)
      } else {
        print("Lists do not match")
        break
      }

      compare$diff <- compare$IntAgeLODDiff - compare$addIntAgeLODDiff
      compare$lod2 <- compare$diff >= 2

      plot <- ggplot(compare, aes(x=IntAgeLODDiff,  y=addIntAgeLODDiff)) +
              geom_point(alpha=0.5) +
              geom_abline(intercept = 0, slope = 1, color="red") +
              geom_abline(intercept = -2, slope = 1, color="blue") +
              scale_x_continuous( name = "LOD score Interactive age pQTL-diff",
                            breaks = seq(0, 15, by = 1),
                            labels = seq(0, 15, by = 1)) +
              scale_y_continuous( name = paste0("LOD score (X | ", med_genename," mRNA)"),
                            breaks = seq(0, 15, by = 1),
                            labels = seq(0, 15, by = 1)) +
              theme_bw() +
              labs(title=paste0("pERK1 upstream mediator: ", med_genename, " (mRNA)"),
                    subtitle = paste0("Chr ", chr, " total: ", nrow(compare), " genes, threshold > 6 \n",
                                  "Genes with LOD drop >= 2: ", nrow(compare[compare$lod2 == TRUE,])))

      # Plot Chr7 LOD scores
      pdf(paste0("./QTLscan/output/plots/pERK1_upstream/mRNA/pERK1_upstream_mRNA_mediator_", med_query,"_", med_genename,"_compare_chr7_thr6.pdf"), width = 9, heigh =9)
      print(plot)
      dev.off()
  }
}
