library(tidyverse)
library(gridExtra)
library(ggsci)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO188b_kidney.RData")
Upheno <- read.delim("./Phenotype/phenotypes/formatted/JAC_CS_urine_chem_v1.txt", sep = "\t")

pheno <- Upheno %>%
         mutate(
           Mouse.ID = mouse.id,
           duplicated = (duplicated(Upheno$mouse.id) | duplicated(Upheno$mouse.id, fromLast = TRUE)),
           Cre = cr.u,
           Alb = ma.u,
           Phs = phs.u) %>%
         filter((Mouse.ID %in% annot.samples$Mouse.ID) & duplicated == FALSE) %>%
         select(Mouse.ID, Alb, Phs, Cre)


annot.samples <- annot.samples[annot.samples$Mouse.ID %in% pheno$Mouse.ID,]
expr.mrna <- expr.mrna[rownames(expr.mrna) %in% pheno$Mouse.ID,]
expr.protein <- expr.protein[rownames(expr.protein) %in% pheno$Mouse.ID,]

# Identify gene name
gene1m <- annot.mrna[annot.mrna$symbol == "Slc34a1",]$id
gene1p <- annot.protein[annot.protein$symbol == "Slc34a1",]$id
gene2m <- annot.mrna[annot.mrna$symbol == "Slc34a3",]$id
gene2p <- annot.protein[annot.protein$symbol == "Slc34a3",]$id
gene3m <- annot.mrna[annot.mrna$symbol == "Lrp2",]$id
gene3p <- annot.protein[annot.protein$symbol == "Lrp2",]$id

df <- annot.samples %>%
      mutate(
        Slc34a1_mRNA = expr.mrna[, gene1m],
        Slc34a1_protein = expr.protein[, gene1p],
        Slc34a3_mRNA = expr.mrna[, gene2m],
        Slc34a3_protein = expr.protein[, gene2p],
        Lrp2_mRNA = expr.mrna[, gene3m],
        Lrp2_protein = expr.protein[, gene3p],
        Sex = as.factor(Sex),
        Generation = factor(Generation, levels = c("G8","G9","G10","G11","G12")),
        Age = factor(as.character(Age), levels = c("6","12","18")),
        Alb = pheno$Alb,
        Phs = pheno$Phs,
        Cre = pheno$Cre) %>%
      select(Mouse.ID, Sex, Generation, Age, Slc34a1_mRNA, Slc34a1_protein, Slc34a3_mRNA, Slc34a3_protein, Lrp2_mRNA, Lrp2_protein, Alb, Phs, Cre)

data <- df %>%
  mutate(AlbCreRatio = Alb/Cre * 100,
        pass30 = AlbCreRatio > 30,
        pass30 = as.factor(pass30)) %>%
  filter(!is.na(pass30))

lrp2_m <- ggplot(data, aes(x= pass30, y = Lrp2_mRNA, colour=Age))+
              geom_boxplot()+
              facet_wrap(~Sex)+
              coord_flip()

lrp2_p <- ggplot(data, aes(x= pass30, y = Lrp2_protein, colour=Age))+
              geom_boxplot()+
              facet_wrap(~Sex)+
              coord_flip()

pdf("ACR_test_plot.pdf", height = 8, width = 10)
grid.arrange(lrp2_m,lrp2_p, ncol =1)
dev.off()



pdf("./Plot/Figure3.pdf", width = 20, height = 10)
grid.arrange(r1c1,r1c2,plot_r1c3,plot_r1c4, r2c1,r2c2,plot_r2c3,plot_r2c4, r3c1,r3c2,plot_r3c3,plot_r3c4, ncol = 4)
dev.off()
