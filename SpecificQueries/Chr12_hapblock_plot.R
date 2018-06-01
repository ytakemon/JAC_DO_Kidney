library(qtl2)
library(qtl2convert)
library(tidyverse)
library(reshape2)
#library(ggsci)
options(dplyr.width = Inf) #override column limit
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")

# get table to use
IntAge_output_file <- "./QTLscan/scanBestMarker_mrna/BestMarker_BestperGene_mrna.csv"
table <- read_csv(IntAge_output_file)
table12 <- table %>% filter(IntAgeChr == 12, IntAgeLODDiff > 7)
PosList <- unique(table12$IntAgePos) %>% sort()
MarkerList <- snps[snps$chr == 12 & snps$bp %in% PosList,]$marker

sub_genoprobs <- genoprobs[,,MarkerList]



testcor <- cor(sub_genoprobs[,,1],sub_genoprobs[,,2])

cor(sub_genoprobs[,,1],sub_genoprobs[,,2])




sub_probs <- probs_doqtl_to_qtl2(sub_genoprobs, snps, pos_column = "bp")

K <- calc_kinship(sub_probs, cores =10)
K <- as.data.frame(melt(K))

pdf("./Plot/Chr12_clustermarker_kinship.pdf", height = 8, width = 8)
qplot(x= Var1, y = Var2, data = K, fill = value, geom="tile") +
scale_fill_gradient(limits=c(0, 1), low ="red", high ="white")
dev.off()
