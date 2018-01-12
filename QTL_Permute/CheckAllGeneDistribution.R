# Check distribution gene expression distibution and pick 2 genes with different
# distributions.
library(ggplot2)
library(dplyr)
library(reshape2)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")

# mrna
identical(colnames(raw.mrna),annot.mrna$id)
mrna <- annot.mrna
mrna$expr <- NA
for (i in 1:nrow(annot.mrna)){
  mrna$expr[i] <- !all(raw.mrna[,i] < 1)
}
expressed <- mrna[mrna$duplicated == FALSE,]
expressed <- expressed[expressed$expr == TRUE,]
mrna <- raw.mrna
mrna <- mrna[,expressed$id]
mrna <- mrna[,expressed[expressed$biotype == "protein_coding",]$id]

sub <- mrna[,c(90:92)]
#See distribution of all genes
df <- melt(sub)
ggplot(df, aes(x = Var2, y = value)) +
  geom_boxplot()

sub <- mrna[,c(90,92)]
#See distribution of all genes
df <- melt(sub)
ggplot(df, aes(x = Var2, y = value)) +
  geom_boxplot()

# check
annot.mrna[annot.mrna$id %in% colnames(mrna)[c(90,91)],]

# Picking:
#> annot.mrna[annot.mrna$id %in% colnames(mrna)[c(90,91)],]
#                    id symbol chr    start      end strand middle_point
#387 ENSMUSG00000002129  Sf3a1  11  4160350  4182541      1      4171446
#390 ENSMUSG00000002204  Napsa   7 44572432 44586862      1     44579647
#    nearest_snp        biotype duplicated
#387       36336 protein_coding      FALSE
#390       24028 protein_coding      FALSE

# protein
