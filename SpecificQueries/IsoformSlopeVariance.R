library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")

load("./RNAseq_data/DO188b_kidney_noprobs.RData")

# Protein: 6716
# Pairs: 6667
# Diff : 49

# annot.mrna: 22312
# expr.mrna : 22243
# diff : 69

# protein
pairs <- annot.protein[annot.protein$gene_id %in% annot.mrna$id,]
nopairs <- annot.protein[!(annot.protein$gene_id %in% annot.mrna$id),]

# What are the proteins that do not have mRNA information?
write.csv(nopairs, file = "./AnnotProt_notin_pair.csv", row.names = FALSE, quote = FALSE)

# What are the slopes of the proteins that are isoforms? do they switch quadrants?
dup <- annot.mrna[annot.mrna$duplicated == TRUE,]
dup_gene <- unique(dup$id)

df <- read.csv("./Anova_output/kidney_anova_slope_output.csv")
df <- select(df, id, gene_id, symbol, m.mRNA_Age.Sex, m.Prot_Age.Sex, p.mRNA_Age.Sex, p.Prot_Age.Sex)

iso <- df[df$gene_id %in% dup_gene,]
iso <- arrange(iso, gene_id)

iso$quadI <- ((iso$m.mRNA_Age.Sex > 0) & (iso$m.Prot_Age.Sex > 0))
iso$quadII <- ((iso$m.mRNA_Age.Sex < 0) & (iso$m.Prot_Age.Sex > 0))
iso$quadIII <- ((iso$m.mRNA_Age.Sex < 0) & (iso$m.Prot_Age.Sex < 0))
iso$quadIV <- ((iso$m.mRNA_Age.Sex > 0) & (iso$m.Prot_Age.Sex < 0))

write.csv(iso, file = "./Isoforms_slopes.csv", row.names = FALSE, quote = FALSE)
