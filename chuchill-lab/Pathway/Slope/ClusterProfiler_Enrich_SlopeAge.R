# R/3.4.4
# http://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

# set up directories
basedir <- "/projects/churchill-lab/projects/JAC/Takemon_DO_crosssectional_kidney/"
setwd(basedir)

# Load data
SlopeResults <- read.csv("./results/ANOVA/Condensed_kidney_anova_output_slope_pvalONLYsig.csv")

# subset for only significant genes < 0.05
sig <- filter(SlopeResults, sig.mRNA_Age.Sex == TRUE & sig.Prot_Age.Sex == TRUE)
# increasing mRNA & protein
SigA <- filter(sig, m.mRNA_Age.Sex > 0 & m.Prot_Age.Sex > 0)
# decreasing mRNA & protein
SigB <- filter(sig, m.mRNA_Age.Sex < 0 & m.Prot_Age.Sex < 0)
# increasing mRNA & decreasing protein
SigC <- filter(sig, m.mRNA_Age.Sex > 0 & m.Prot_Age.Sex < 0)
# decreasing mRNA & increasing protein
SigD <- filter(sig, m.mRNA_Age.Sex < 0 & m.Prot_Age.Sex > 0)


# Get universe list:
universe_df <- bitr(as.character(SlopeResults$gene_id),
                    fromType = "ENSEMBL",
                    toType = c("ENTREZID", "SYMBOL"),
                    OrgDb = org.Mm.eg.db)

# Warning message:
# In bitr(as.character(SlopeResults$gene_id), fromType = "ENSEMBL",  :
#  1.83% of input gene IDs are fail to map...

QuadList <- list(SigA, SigB, SigC, SigD)
for(i in 1:length(QuadList)){
  Sigdf <- QuadList[[i]]
  gene_df <- bitr(as.character(Sigdf$gene_id),
                    fromType = "ENSEMBL",
                    toType = c("ENTREZID", "SYMBOL"),
                    OrgDb = org.Mm.eg.db)

  # GO over-representation test
  # GO:bp
  go_bp <- enrichGO( gene = gene_df$ENTREZID,
                   universe = universe_df$ENTREZID,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

  # GO:cc
  go_cc <- enrichGO( gene = gene_df$ENTREZID,
                   universe = universe_df$ENTREZID,
                   OrgDb = org.Mm.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)

  # KEGG over-representation test
  kegg <- enrichKEGG( gene = gene_df$ENTREZID,
                        organism = "mmu",
                        pvalueCutoff = 0.05)

  # Save enrich KEGG, but first convert all ENTREZID into symbol for readability.
  if(!nrow(kegg) == 0){
    kegg_result <- kegg@result
    kegg_result$geneID_symbol <- NA
    for(n in 1:nrow(kegg_result)) {
      print(n)
      list <- kegg_result$geneID[n]
      list <- strsplit(list, "/")[[1]]
      list <- as.data.frame(list)
      list$ENSEMBL <- NA
      for(s in 1:nrow(list)){
        list$ENSEMBL[s] <- universe_df[universe_df$ENTREZID %in% list[s,"list"],]$SYMBOL[1]
      }

      kegg_result$geneID_symbol[n] <- paste(list$ENSEMBL, collapse = "/")
    }
  }else{
    kegg_result <- kegg@result
  }


  # Save files
  if(i == 1){
    type <- "IncRNA_IncProt"
  } else if (i == 2){
    type <- "DecRNA_DecProt"
  } else if (i == 3){
    type <- "IncRNA_DecProt"
  } else if (i == 4){
    type <- "DecRNA_IncProt"
  }

  go_result <- go_bp@result
  write.csv(go_result, paste0("./results/Pathways/AgeSlope_EnrichGObp_",type,".csv"), row.names=FALSE)

  go_result <- go_cc@result
  write.csv(go_result, paste0("./results/Pathways/AgeSlope_EnrichGOcc_",type,".csv"), row.names=FALSE)

  write.csv(kegg_result, paste0("./results/Pathways/AgeSlope_EnrichKEGG_",type,".csv"), row.names=FALSE)
}
