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

# check protein expression
select <- annot.protein[annot.protein$gene_id %in% colnames(mrna)[c(90,91)],]
prot <- raw.protein
sub <- prot[,select$id]
df <- melt(sub)
ggplot(df, aes(Var2, y = value)) +
  geom_boxplot()

# Doesn't look very different. Picking different proteins
identical(colnames(raw.protein),annot.protein$id)
expressed <- annot.protein[annot.protein$biotype == "protein_coding",]
prot <- raw.protein[, annot.protein$id]

sub <- prot[,146:148]
df <- melt(sub)
ggplot(df, aes(x = Var2, y = value)) +
  geom_boxplot()

# check
annot.protein[annot.protein$id %in% colnames(prot)[c(146,148)],]

# Picking for protein:
#> annot.protein[annot.protein$id %in% colnames(prot)[c(146,148)],]
#                    id            gene_id symbol chr    start      end strand
#197 ENSMUSP00000003529 ENSMUSG00000003437   Paf1   7 28392951 28399388      1
#200 ENSMUSP00000003569 ENSMUSG00000003477   Inmt   6 55170626 55175043     -1
#    middle_point nearest_snp        biotype
#197     28396170       23520 protein_coding
#200     55172834       20641 protein_coding
