# R/3.4.1
library(dplyr)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE) # args <- "kidney_anova_slope_output.csv"
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("RNAseq_data/DO188b_kidney_noprobs.RData")
output <- list.files(path = "./Anova_output/", pattern = paste0("^",args[[1]]), recursive = TRUE)
data <- read.csv(paste0("./Anova_output/",output[[1]]), header = T)

# For only slopes and sigcols
mcols <- grep("^m.", colnames(data))
sigcols <- grep("^p.", colnames(data))
data <- data[, c(1:8, mcols, sigcols)]

# Age (not interacting with sex)------------------------------------------------
# Using only age realated mRNA and Protein that are significant < 0.05
df <- data[data$p.mRNA_Age.Sex <= 0.05,]
df <- df[df$p.Prot_Age.Sex <= 0.05,]
# N = 693

# Count number of genes in each quadrant
# mRNA-Protein
total <- nrow(df)
quadI <- nrow(df[((df$m.mRNA_Age.Sex > 0) & (df$m.Prot_Age.Sex > 0)),])
quadIV <- nrow(df[((df$m.mRNA_Age.Sex > 0) & (df$m.Prot_Age.Sex < 0)),])
quadII <- nrow(df[((df$m.mRNA_Age.Sex < 0) & (df$m.Prot_Age.Sex > 0)),])
quadIII <- nrow(df[((df$m.mRNA_Age.Sex < 0) & (df$m.Prot_Age.Sex < 0)),])
# Export gene list in each quadrant
Gene_list <- df
Gene_list$quadI <- ((Gene_list$m.mRNA_Age.Sex > 0) & (Gene_list$m.Prot_Age.Sex > 0))
Gene_list$quadII <- ((df$m.mRNA_Age.Sex < 0) & (df$m.Prot_Age.Sex > 0))
Gene_list$quadIII <- ((df$m.mRNA_Age.Sex < 0) & (df$m.Prot_Age.Sex < 0))
Gene_list$quadIV <- ((df$m.mRNA_Age.Sex > 0) & (df$m.Prot_Age.Sex < 0))
write.csv(Gene_list, "./Anova_output/gene_lists/Quad_Age.csv", row.names = FALSE, quote = FALSE)

# Plot age realted mRNA/Protein slopes
pdf("./Plot/slope_mRNA_Prot_Age.pdf", width = 6, heigh =6)
ggplot(df, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex)) +
      geom_point() +
      scale_x_continuous("mRNA-Age.Sex",
                         breaks = seq(-1, 1, by = 0.25)) +
      scale_y_continuous("Protein-Age.Sex",
                         breaks = seq(-1, 1, by = 0.25)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs( title = "Direction of change in gene expression with Age",
            subtitle = paste0("Slope: RNA v. Protein, N = ", total)) +
      annotate("text", x = 0.25, y = 0.85, label = paste("I, N = ", quadI)) +
      annotate("text", x = 0.25, y = -0.85, label = paste("IV, N = ", quadIV)) +
      annotate("text", x = -0.25, y = 0.85, label = paste("II, N = ", quadII)) +
      annotate("text", x = -0.25, y = -0.85, label = paste("III, N = ", quadIII)) +
      theme(panel.background = element_blank(),
            panel.border = element_rect( colour = "black", fill = NA))
dev.off()

# Age interacting with Sex -----------------------------------------------------
# Using only age realated mRNA and Protein that are significant < 0.05
df <- data[data$p.mRNA_Interaction <= 0.05,]
df <- df[df$p.Prot_Interaction <= 0.05,]
# N = 168

# Count number of genes in each quadrant
# mRNA-Protein
total <- nrow(df)
quadI <- nrow(df[((df$m.mRNA_Interaction > 0) & (df$m.Prot_Interaction > 0)),])
quadIV <- nrow(df[((df$m.mRNA_Interaction > 0) & (df$m.Prot_Interaction < 0)),])
quadII <- nrow(df[((df$m.mRNA_Interaction < 0) & (df$m.Prot_Interaction > 0)),])
quadIII <- nrow(df[((df$m.mRNA_Interaction < 0) & (df$m.Prot_Interaction < 0)),])
# Export gene list in each quadrant
Gene_list <- df
Gene_list$quadI <- ((Gene_list$m.mRNA_Interaction > 0) & (Gene_list$m.Prot_Interaction > 0))
Gene_list$quadII <- ((df$m.mRNA_Interaction < 0) & (df$m.Prot_Interaction > 0))
Gene_list$quadIII <- ((df$m.mRNA_Interaction < 0) & (df$m.Prot_Interaction < 0))
Gene_list$quadIV <- ((df$m.mRNA_Interaction > 0) & (df$m.Prot_Interaction < 0))
write.csv(Gene_list, "./Anova_output/gene_lists/Quad_AgeIntSex.csv", row.names = FALSE, quote = FALSE)
# Plot age realted mRNA/Protein slopes
pdf("./Plot/slope_mRNA_Prot_AgeIntSex.pdf", width = 6, heigh =6)
ggplot(df, aes(x = m.mRNA_Interaction, y = m.Prot_Interaction)) +
      geom_point() +
      scale_x_continuous("mRNA-Age.Sex",
                         breaks = seq(-1, 1, by = 0.02)) +
      scale_y_continuous("Protein-Age.Sex",
                         breaks = seq(-1, 1, by = 0.02)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs( title = "Direction of change in gene expression with Age*Sex",
            subtitle = paste0("Slope: RNA v. Protein, N = ", total)) +
      annotate("text", x = 0.06, y = 0.14, label = paste("I, N = ", quadI)) +
      annotate("text", x = 0.06, y = -0.14, label = paste("IV, N = ", quadIV)) +
      annotate("text", x = -0.06, y = 0.14, label = paste("II, N = ", quadII)) +
      annotate("text", x = -0.06, y = -0.14, label = paste("III, N = ", quadIII)) +
      theme(panel.background = element_blank(),
            panel.border = element_rect( colour = "black", fill = NA))
dev.off()
