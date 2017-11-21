library(ggplot2)
library(dplyr)
library(gridExtra)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")
erk1 <- read.delim("./Phenotype/phenotypes/JAC_WB_kidney_ERK.txt")
pheno <- "Phospho_ERK1_ratio"

# Cleanup data and subset to match samples
erk1 <- erk1[complete.cases(erk1[,pheno]),]
erk1 <- erk1[erk1$ID %in% annot.samples$Mouse.ID, ]
erk1$duplicated <- (duplicated(erk1$ID) | duplicated(erk1$ID, fromLast = TRUE))
erk1 <- erk1[erk1$duplicated == FALSE,]
erk1 <- erk1[erk1$ID %in% rownames(annot.samples),] # 55 animals (33 females, 22 males)
rownames(erk1) <- erk1$ID

# new df
df <- annot.samples[,1:4]
df <- df[rownames(df) %in% rownames(erk1),]
df$pheno <- log(erk1[,pheno])
df$Age <- as.factor(df$Age)
df$Sex <- as.factor(df$Sex)

fitF <- df[df$Sex == "F",]
fitF <- lm(pheno ~ Age, fitF)
int <- signif(coef(fitF)[1], 3)
slope <- signif(coef(fitF)[2], 3)
pval <- signif(summary(fitF)$coefficients[2,4], 3)
r2 <- signif(summary(fitF)$adj.r.squared, 3)
eqF <- paste0("f(female) = ", slope, "x", int, ", R^2 = ", r2, ", pval = ", pval)

fitM <- df[df$Sex == "M",]
fitM <- lm(pheno ~ Age, fitM)
int <- signif(coef(fitM)[1], 3)
slope <- signif(coef(fitM)[2], 3)
pval <- signif(summary(fitM)$coefficients[2,4], 3)
r2 <- signif(summary(fitM)$adj.r.squared, 3)
eqM <- paste0("f(male) = ", slope, "x", int, ", R^2 = ", r2, ", pval = ", pval)

erk1 <- ggplot(df, aes(x = Age, y = pheno, fill = Sex))+
              geom_smooth(method = "lm", se=FALSE, aes(group = Sex, colour = Sex, fill = Sex)) +
              geom_point(aes(colour = Sex, fill = Sex, alpha = 0.5), size = 2, shape = 21, position = position_jitterdodge())+
              labs( title = "p-ERK1 / ERK1 ratio",
                    subtitle = paste(eqF, "\n", eqM),
                    y = "log-transformed p-ERK1 / ERK1 ratio",
                    x = "Age") +
              guides(alpha = FALSE) +
              theme(panel.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(color = "black"))

pdf(paste0("./Plot/pERK_ratio_Age.pdf"), width = 8, height = 6)
erk1
dev.off()
