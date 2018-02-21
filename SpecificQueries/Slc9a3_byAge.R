library(dplyr)
library(ggsci)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")

Slc_df <- annot.samples
Slc_df$Age <- as.factor(as.character(Slc_df$Age))
Slc_df$Age <- factor(Slc_df$Age, levels(Slc_df$Age)[c(3,1,2)]) # Change levels order


# find gene id
id <- annot.mrna[annot.mrna$symbol == "Slc9a3",]$id
pid <- annot.protein[annot.protein$symbol == "Slc9a3",]$id
Slc_df$Slc9a3_mrna <- expr.mrna[, id]
Slc_df$Slc9a3_mrna_raw <- raw.mrna[, id]
Slc_df$Slc9a3_prot <- expr.protein[, pid]
Slc_df$Slc9a3_prot_raw <- raw.protein[, pid]


pdf("./Plot/Slc9a3cor_byAge.pdf", width = 12, height = 7)
ggplot(Slc_df, aes(x = Slc9a3_mrna_raw, y = Slc9a3_prot_raw, colour = Age)) +
      geom_smooth(method = "lm", se = FALSE, aes(group = Age, colour = Age)) +
      geom_point() +
      theme_bw() +
      labs( title = "Slc9a3 mRNA v. Protein Expression by Age",
            subtitle = paste0("Female: 93 (6mo=33, 12mo=31, 18mo=29) \n Male: 95(6mo=30, 12mo=31, 18mo=34)"),
            y = "SLC9A3 expression (protein)",
            x = "Slc9a3 expression (mRNA)") +
      facet_grid(. ~ Sex) +
      scale_colour_aaas()
dev.off()
