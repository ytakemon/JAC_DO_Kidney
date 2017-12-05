library(dplyr)
library(gridExtra)

setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")
erk1 <- read.delim("./Phenotype/phenotypes/JAC_WB_kidney_ERK.txt", stringsAsFactors = FALSE)
pheno <- "Total_ERK1"

prot <- "Mapk3"
prot <- annot.protein[annot.protein$symbol == prot,]

erk1 <- erk1[erk1$ID %in% rownames(annot.samples),]
annot.samples <- annot.samples[rownames(annot.samples) %in% erk1$ID, ]
expr.protein <- expr.protein[rownames(expr.protein) %in% erk1$ID, ]

identical(erk1$ID, rownames(annot.samples))
identical(erk1$ID, rownames(expr.protein))

erk1$Mapk3 <- expr.protein[, prot$id]

fit <- lm(log(Total_ERK1) ~ Mapk3, erk1)
int <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(summary(fit)$coefficients[2,4], 3)
r2 <- signif(summary(fit)$adj.r.squared, 3)
eq <- paste0("f() = ", slope, "x + ", int, ", R^2 = ", r2, ", pval = ", pval)
plot1 <- ggplot(erk1, aes(x = Mapk3, y = log(Total_ERK1))) +
      geom_point() +
      geom_smooth( method = "lm", se = FALSE) +
      theme_bw() +
      labs(title = "MAPK3 (proteome) vs log(Total_ERK1)",
           subtitle = eq)

fit <- lm(log(Phospho_ERK1) ~ Mapk3, erk1)
int <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(summary(fit)$coefficients[2,4], 3)
r2 <- signif(summary(fit)$adj.r.squared, 3)
eq <- paste0("f() = ", slope, "x + ", int, ", R^2 = ", r2, ", pval = ", pval)
plot2 <- ggplot(erk1, aes(x = Mapk3, y = log(Phospho_ERK1))) +
      geom_point() +
      geom_smooth( method = "lm", se = FALSE) +
      theme_bw() +
      labs(title = "MAPK3 (proteome) vs log(Phospho_ERK1)",
           subtitle = eq)

fit <- lm(log(Total_ERK1) ~ log(Phospho_ERK1), erk1)
int <- signif(coef(fit)[1], 3)
slope <- signif(coef(fit)[2], 3)
pval <- signif(summary(fit)$coefficients[2,4], 3)
r2 <- signif(summary(fit)$adj.r.squared, 3)
eq <- paste0("f() = ", slope, "x + ", int, ", R^2 = ", r2, ", pval = ", pval)
plot3 <- ggplot(erk1, aes(x =log(Phospho_ERK1), y = log(Total_ERK1))) +
      geom_point() +
      geom_smooth( method = "lm", se = FALSE) +
      theme_bw() +
      labs(title = " log(Total_ERK1) vs log(Phospho_ERK1)",
           subtitle = eq)

pdf(paste0("./Plot/Erk_comparisons.pdf"), width = 12, height = 6)
grid.arrange(plot1, plot2, plot3, ncol = 3)
dev.off()
