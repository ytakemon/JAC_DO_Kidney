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
df_age <- data[data$p.mRNA_Age.Sex <= 0.05,]
df_age <- df_age[df_age$p.Prot_Age.Sex <= 0.05,]

# Change in mRNA while protein stays constant ----------------------------------
cols <- c("m.mRNA_Age.Sex", "m.Prot_Age.Sex", "m.mRNA_Interaction", "m.Prot_Interaction",
          "p.mRNA_Age.Sex", "p.Prot_Age.Sex", "p.mRNA_Interaction", "p.Prot_Interaction")
df_deltaAge <- data[data$p.mRNA_Age.Sex <= 0.05,]
df_deltaAge <- df_deltaAge[df_deltaAge$p.Prot_Age.Sex > 0.05,]
sub <- which(colnames(df_deltaAge) %in% cols)
df_deltaAge <- df_deltaAge[, c(1:8,sub)]

# Change in protein while mRNA stays constant ----------------------------------
df_deltaProt_inc <- data[data$p.Prot_Age.Sex <= 1e-10,]
df_deltaProt_inc <- df_deltaProt_inc[df_deltaProt_inc$p.mRNA_Age.Sex > 0.05,]

sub <- which(colnames(df_deltaProt_inc) %in% cols)
df_deltaProt_inc <- df_deltaProt_inc[, c(1:8,sub)]


# Plot total with subset overlay ----------------------------------------------
# Age
total <- nrow(data)
pdf("./Plot/Total_slope_mRNA_Prot_Age.pdf", width = 6, heigh =6)
ggplot() +
      geom_point(data = data, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex, alpha = 0.2), color = "grey") +
      geom_point(data = df_age, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex), colour = "#FF0000") +
      geom_point(data = df_deltaAge, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex), colour = "#06F030") +
      geom_point(data = df_deltaProt_inc, aes(x = m.mRNA_Age.Sex, y = m.Prot_Age.Sex), colour = "#06F030") +
      scale_x_continuous("mRNA-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      scale_y_continuous("Protein-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs( title = "Direction of change in gene expression with Age",
            subtitle = paste0("Slope: mRNA v. Protein, Total = ", total)) +
      guides(colour = FALSE, alpha = FALSE) +
      theme(panel.background = element_blank(),
            panel.border = element_rect( colour = "black", fill = NA))
dev.off()
