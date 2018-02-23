library(ggplot2)
library(dplyr)
library(gridExtra)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO188b_kidney.RData")
Upheno <- read.delim("./Phenotype/phenotypes/formatted/JAC_CS_urine_chem_v1.txt", sep = "\t")

pheno <- Upheno %>%
         mutate(
           Mouse.ID = mouse.id,
           duplicated = (duplicated(Upheno$mouse.id) | duplicated(Upheno$mouse.id, fromLast = TRUE)),
           Alb = ma.u,
           Phs = phs.u) %>%
         filter((Mouse.ID %in% annot.samples$Mouse.ID) & duplicated == FALSE) %>%
         select(Mouse.ID, Alb, Phs)


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
        Phs = pheno$Phs) %>%
      select(Mouse.ID, Sex, Generation, Age, Slc34a1_mRNA, Slc34a1_protein, Slc34a3_mRNA, Slc34a3_protein, Lrp2_mRNA, Lrp2_protein, Alb, Phs)

r1c1 <- ggplot(df, aes(Age, y = Slc34a1_mRNA, colour = Sex)) +
          geom_smooth(method = "lm", se = FALSE, aes(group = Sex, colour = Sex)) +
          geom_point(position = position_jitterdodge(0.4)) +
          theme_bw() +
          labs(y = "Slc34a1 mRNA",
               x = "Age") +
          theme(panel.background = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(color = "black")) +
          scale_colour_aaas()

r1c2 <- ggplot(df, aes(Age, y = Slc34a1_protein, colour = Sex)) +
          geom_smooth(method = "lm", se = FALSE, aes(group = Sex, colour = Sex)) +
          geom_point(position = position_jitterdodge(0.4)) +
          theme_bw() +
          labs(y = "Slc34a1 protein",
               x = "Age") +
          theme(panel.background = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(color = "black")) +
          scale_colour_aaas()

r1c3 <- ggplot(df, aes(x = Slc34a1_mRNA, y = Slc34a1_protein, colour = Age)) +
          geom_smooth(method = "lm", se = FALSE, aes(group = Age, colour = Age)) +
          geom_point() +
          theme_bw() +
          labs(y = "Slc34a1 protein",
               x = "Slc34a1 mRNA") +
          theme(panel.background = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                axis.line = element_line(color = "black")) +
          facet_grid(. ~ Sex) +
          scale_colour_aaas()

r1c4 <- ggplot(df, aes(x = Slc34a1_protein, y = Phs, colour = Age)) +
          geom_smooth(method = "lm", se = FALSE, aes(group = Age, colour = Age)) +
          geom_point() +
          coord_cartesian(ylim = c(0,1000))+
          theme_bw() +
          labs(y = "Phsosphate",
               x = "Slc34a1 mRNA") +
          theme(panel.background = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                axis.line = element_line(color = "black")) +
          facet_grid(. ~ Sex) +
          scale_colour_aaas()

r2c1 <- ggplot(df, aes(Age, y = Slc34a3_mRNA, colour = Sex)) +
          geom_smooth(method = "lm", se = FALSE, aes(group = Sex, colour = Sex)) +
          geom_point(position = position_jitterdodge(0.4)) +
          theme_bw() +
          labs(y = "Slc34a3 mRNA",
               x = "Age") +
          theme(panel.background = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(color = "black")) +
          scale_colour_aaas()

r2c2 <- ggplot(df, aes(Age, y = Slc34a3_protein, colour = Sex)) +
          geom_smooth(method = "lm", se = FALSE, aes(group = Sex, colour = Sex)) +
          geom_point(position = position_jitterdodge(0.4)) +
          theme_bw() +
          labs(y = "Slc34a3 protein",
               x = "Age") +
          theme(panel.background = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(color = "black")) +
          scale_colour_aaas()

r2c3 <- ggplot(df, aes(x = Slc34a3_mRNA, y = Slc34a3_protein, colour = Age)) +
          geom_smooth(method = "lm", se = FALSE, aes(group = Age, colour = Age)) +
          geom_point() +
          theme_bw() +
          labs(y = "Slc34a3 protein",
               x = "Slc34a3 mRNA") +
          theme(panel.background = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                axis.line = element_line(color = "black")) +
          facet_grid(. ~ Sex) +
          scale_colour_aaas()

r2c4 <- ggplot(df, aes(x = Slc34a3_protein, y = Phs, colour = Age)) +
          geom_smooth(method = "lm", se = FALSE, aes(group = Age, colour = Age)) +
          geom_point() +
          coord_cartesian(ylim = c(0,1000))+
          theme_bw() +
          labs(y = "Phosphate",
               x = "Slc34a3 mRNA") +
          theme(panel.background = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                axis.line = element_line(color = "black")) +
          facet_grid(. ~ Sex) +
          scale_colour_aaas()

r3c1 <- ggplot(df, aes(Age, y = Lrp2_mRNA, colour = Sex)) +
          geom_smooth(method = "lm", se = FALSE, aes(group = Sex, colour = Sex)) +
          geom_point(position = position_jitterdodge(0.4)) +
          theme_bw() +
          labs(y = "Lrp2 mRNA",
               x = "Age") +
          theme(panel.background = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(color = "black")) +
          scale_colour_aaas()

r3c2 <- ggplot(df, aes(Age, y = Lrp2_protein, colour = Sex)) +
          geom_smooth(method = "lm", se = FALSE, aes(group = Sex, colour = Sex)) +
          geom_point(position = position_jitterdodge(0.4)) +
          theme_bw() +
          labs(y = "Lrp2 protein",
               x = "Age") +
          theme(panel.background = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                panel.border = element_blank(),
                axis.line = element_line(color = "black")) +
          scale_colour_aaas()

r3c3 <- ggplot(df, aes(x = Lrp2_mRNA, y = Lrp2_protein, colour = Age)) +
          geom_smooth(method = "lm", se = FALSE, aes(group = Age, colour = Age)) +
          geom_point() +
          theme_bw() +
          labs(y = "Lrp2 protein",
               x = "Lrp2 mRNA") +
          theme(panel.background = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                axis.line = element_line(color = "black")) +
          facet_grid(. ~ Sex) +
          scale_colour_aaas()

r3c4 <- ggplot(df_new, aes(x = Lrp2_protein, y = Alb, colour = Age)) +
          geom_smooth(method = "lm", se = FALSE, aes(group = Age, colour = Age)) +
          geom_point() +
          scale_y_continuous(limits = c(0,8))+
          theme_bw() +
          labs(y = "Albumin",
               x = "Lrp2 protein") +
          theme(panel.background = element_blank(),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                axis.line = element_line(color = "black")) +
          facet_grid(. ~ Sex) +
          scale_colour_aaas()


pdf("./Plot/Figure3.pdf", width = 20, height = 10)
grid.arrange(r1c1,r1c2,r1c3,r1c4, r2c1,r2c2,r2c3,r2c4, r3c1,r3c2,r3c3,r3c4, ncol = 4)
dev.off()
