# R/3.4.4
# http://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

# set up directories
basedir <- "/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/"
setwd(basedir)

# Load data
load("./RNAseq_data/DO188b_kidney.RData")
rm(genoprobs, G, Glist, N, raw.mrna, snps)
AnovResults <- read.csv("./Anova_output/protein.kidney_anova_table.csv")

# subset for only significant genes < 0.05
sig <- AnovResults[AnovResults$sig.Prot_Age.Sex == TRUE,]

# GO over-representation test needs ENTREZID as input
gene_df <- bitr(as.character(sig$gene_id),
                fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Mm.eg.db)
# Somthing to consider:
# Warning message:
# In bitr(as.character(sig$gene_id), fromType = "ENSEMBL", toType = c("ENTREZID",  :
#  1.7% of input gene IDs are fail to map...

universe_df <- bitr(as.character(AnovResults$gene_id),
                 fromType = "ENSEMBL",
                 toType = c("ENTREZID", "SYMBOL"),
                 OrgDb = org.Mm.eg.db)
# Somthing to consider:
# Warning message:
# In bitr(as.character(AnovResults$gene_id), fromType = "ENSEMBL",  :
#  1.85% of input gene IDs are fail to map...

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

go_result <- go_bp@result
write.csv(go_result, "./Pathways/PROT_EnrichGObp_Age.csv", row.names = FALSE)

# GO:cc
go_cc <- enrichGO( gene = gene_df$ENTREZID,
                 universe = universe_df$ENTREZID,
                 OrgDb = org.Mm.eg.db,
                 ont = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05,
                 readable = TRUE)

go_result <- go_cc@result
write.csv(go_result, "./Pathways/PROT_EnrichGOcc_Age.csv", row.names = FALSE)

# KEGG over-representation test
kegg <- enrichKEGG( gene = gene_df$ENTREZID,
                      organism = "mmu",
                      pvalueCutoff = 0.05)

# Save enrich KEGG, but first convert all ENTREZID into symbol for readability.
kegg_result <- kegg@result
kegg_result$geneID_symbol <- NA
for(i in 1:nrow(kegg_result)) {
  print(i)
  list <- kegg_result$geneID[i]
  list <- strsplit(list, "/")[[1]]
  list <- as.data.frame(list)
  list$ENSEMBL <- NA
  for(n in 1:nrow(list)){
    list$ENSEMBL[n] <- universe_df[universe_df$ENTREZID %in% list[n,"list"],]$SYMBOL[1]
  }

  kegg_result$geneID_symbol[i] <- paste(list$ENSEMBL, collapse = "/")
}
write.csv(kegg_result, "./Pathways/PROT_EnrichKEGG_Age.csv", row.names = FALSE)
