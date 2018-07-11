# R/3.4.4
# http://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

# set up directories
basedir <- "/projects/korstanje-lab/ytakemon/JAC_DO_Kgene_idney/"
setwd(basedir)

# Define interst:
type <- "protein"
chr <- "7"
thr <- 7

# Load data
if(type == "mrna"){
  ScanResult <- read.csv("./QTLscan/scanBestMarker_mrna/BestMarker_BestperGene_mrna.csv", stringsAsFactors = FALSE)
} else if(type == "protein"){
  ScanResult <- read.csv("./QTLscan/scanBestMarker_protein/BestMarker_BestperGene_protein.csv", stringsAsFactors = FALSE)
}

# subset for only genes in chromosome:
chrList <- ScanResult[ScanResult$IntAgeChr == chr & ScanResult$IntAgeLODDiff > thr,]

# Get universe list:
universe_df <- bitr(as.character(ScanResult$gene_id),
                    fromType = "ENSEMBL",
                    toType = c("ENTREZID", "SYMBOL"),
                    OrgDb = org.Mm.eg.db)

# Warning message:
# In bitr(as.character(ScanResult$gene_id), fromType = "ENSEMBL",  :
#  1.83% of input gene IDs are fail to map...

  gene_df <- bitr(as.character(chrList$gene_id),
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

  go_result <- go_bp@result
  write.csv(go_result, paste0("./Pathways/",type,"QTLscan_chr",chr,"_thr",thr,"_EnrichGObp.csv"), row.names=FALSE)

  go_result <- go_cc@result
  write.csv(go_result, paste0("./Pathways/",type,"QTLscan_chr",chr,"_thr",thr,"_EnrichGOcc.csv"), row.names=FALSE)

  write.csv(kegg_result, paste0("./Pathways/",type,"QTLscan_chr",chr,"_thr",thr,"_EnrichKEGG.csv"), row.names=FALSE)
