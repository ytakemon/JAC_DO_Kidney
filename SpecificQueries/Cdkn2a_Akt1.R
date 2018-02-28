library(ggplot2)
library(dplyr)
library(ggsci)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")


df <- annot.samples
df$Age <- as.factor(as.character(df$Age))
df$Age <- factor(df$Age, levels(df$Age)[c(3,1,2)]) # Change levels order

# find gene id
id <- annot.mrna[annot.mrna$symbol == "Cdkn2a",]$id
aktid <- annot.mrna[annot.mrna$symbol == "Akt1",]$id
df$Cdkn2a_mrna <- expr.mrna[, id]
df$Akt1_mrna <- expr.mrna[, aktid]

pdf("./Plot/Cdkn2a_Akt1_byAge.pdf", width = 12, height = 7)
ggplot(df, aes(x = Cdkn2a_mrna, y = Akt1_mrna, colour = Age)) +
      geom_smooth(method = "lm", se = FALSE, aes(group = Age, colour = Age)) +
      geom_point() +
      theme_bw() +
      labs( title = "Cdkn2a mRNA v. Akt1 mRNA by Age",
            subtitle = paste0("Female: 93 (6mo=33, 12mo=31, 18mo=29) \n Male: 95(6mo=30, 12mo=31, 18mo=34)"),
            y = "Akt1 expression (Rank normalized mRNA)",
            x = "Cdkn2a expression (Rank normalized mRNA)") +
      facet_grid(. ~ Sex) +
      scale_colour_aaas()
dev.off()

erk1 <- read.delim("./Phenotype/phenotypes/JAC_WB_kidney_ERK.txt")
erk1 <- erk1 %>% filter(erk1$ID %in% df$Mouse.ID & !is.na(erk1$Phospho_ERK1_ratio))
df <- df %>% filter(df$Mouse.ID %in% erk1$ID)
df$pERKratio <- erk1$Phospho_ERK1_ratio
count <- table(df$Age, df$Sex)

pdf("./Plot/Cdkn2a_pErk1ratio_byAge.pdf", width = 12, height = 7)
ggplot(df, aes(x = Cdkn2a_mrna, y = log(pERKratio), colour = Age)) +
      geom_smooth(method = "lm", se = FALSE, aes(group = Age, colour = Age)) +
      geom_point() +
      theme_bw() +
      labs( title = "Cdkn2a mRNA v. pERK1 mRNA by Age",
            subtitle = paste0("Female: ",sum(count[,"F"])," (6mo=",count["6","F"],", 12mo=",count["12","F"],", 18mo=",count["18","F"],") \n Male: ",sum(count[,"M"]),"(6mo=",count["6","M"],", 12mo=",count["12","M"],", 18mo=",count["18","M"],")"),
            y = "pERK1 expression (log normalized)",
            x = "Cdkn2a expression (Rank normalized mRNA)") +
      facet_grid(. ~ Sex) +
      scale_colour_aaas()
dev.off()
