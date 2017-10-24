library(ggplot2)
library(dplyr)
library(broom)
library(GGally)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO188b_kidney.RData")
pheno <- read.delim("./Phenotype/phenotypes/formatted/JAC_CS_urine_chem_v1.txt", sep = "\t")

# Subset pheno to match count for animal id
pheno$duplicated <- (duplicated(pheno$mouse.id) | duplicated(pheno$mouse.id, fromLast = TRUE))
pheno <- pheno[pheno$duplicated == FALSE,]
pheno <- pheno[pheno$mouse.id %in% annot.samples$Mouse.ID,] # only 176 animals with phenotypes
pheno <- arrange(pheno, mouse.id)
pheno$phs.cr.u <- (pheno$phs.u / pheno$cr.u) * 100
annot.samples <- annot.samples[annot.samples$Mouse.ID %in% pheno$mouse.id,]
expr.mrna <- expr.mrna[rownames(expr.mrna) %in% pheno$mouse.id,]
expr.protein <- expr.protein[rownames(expr.protein) %in% pheno$mouse.id,]

# Identify gene name
gene1 <- "Slc34a1"
gene2 <- "Slc34a3"
other.ids <- function(gene.name, level) {
  if (level == "mRNA") {
    sel <- which(mRNA.list$symbol == gene.name)[1]
    if (!is.na(sel)) return(mRNA.list[sel,]) else return(c(NA,NA,NA))
  }
  if (level == "protein") {
    sel <- which(protein.list$symbol == gene.name)[1]
    if (!is.na(sel)) return(protein.list[sel,]) else return(c(NA,NA,NA))
  }
}
gene1 <- other.ids(gene1, "mRNA")
gene2 <- other.ids(gene2, "mRNA")

# new df
df <- annot.samples[,1:4]
df <- cbind(df, expr.mrna[, gene1$id], expr.mrna[, gene2$id], expr.protein[, gene1$protein_id], expr.protein[, gene2$protein_id], pheno$phs.cr.u)
colnames(df)[5:9] <- c("Slc34a1_mRNA", "Scl34a3_mRNA", "Slc34a1_protein", "Slc34a3_protein", "Phs.Cre.U")
df$Age <- as.factor(df$Age)


upper_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) +
       theme(panel.background = element_blank(),
             panel.grid.minor = element_blank(),
             panel.grid.major = element_blank(),
             panel.border = element_blank(),
             axis.line = element_line(color = "black"))
  p
}
diag_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) +
       geom_density()+
       theme(panel.background = element_blank(),
             panel.grid.minor = element_blank(),
             panel.grid.major = element_blank(),
             panel.border = element_blank(),
             axis.line = element_line(color = "black"))
  p
}
lower_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) +
       geom_point() +
       geom_smooth(method = lm, se = FALSE, ...) +
       theme(panel.background = element_blank(),
             panel.grid.minor = element_blank(),
             panel.grid.major = element_blank(),
             panel.border = element_blank(),
             axis.line = element_line(color = "black"))
  p
}

pdf("./Plot/Slc34a13_PhsCre.pdf", height = 7, width = 7)
ggpairs(df,
        mapping = aes(color = Age, alpha = 0.2),
        columns = c("Slc34a1_mRNA", "Scl34a3_mRNA", "Slc34a1_protein", "Slc34a3_protein", "Phs.Cre.U"),
        upper = list( discrete = upper_fn),
        diag = list (continuous = diag_fn),
        lower = list(continuous = lower_fn))
dev.off()

# Test Pearson cor for p-values
# Total
pval <- data.frame( Slc34a1_mRNA = integer(),
                    Scl34a3_mRNA = integer(),
                    Slc34a1_protein = integer(),
                    Slc34a3_protein = integer(),
                    Phs.Cre.U = integer())
col <- c("Slc34a1_mRNA", "Scl34a3_mRNA", "Slc34a1_protein", "Slc34a3_protein", "Phs.Cre.U")

for (r in 1:length(col)){
  for (c in 1:length(col)){
      pval[r,c] <- cor.test(df[,(4+r)], df[,(4+c)], alternative = "two.sided", method = "pearson")$p.value
  }
}
rownames(pval) <- col
write.csv(pval, file = "./SpecificQ/Slc34a13_PhsCre_pval.csv", row.names = FALSE, quote = FALSE)

# 6 mo.
df6 <- df[df$Age == 6,]
pval <- data.frame( Slc34a1_mRNA = integer(),
                    Scl34a3_mRNA = integer(),
                    Slc34a1_protein = integer(),
                    Slc34a3_protein = integer(),
                    Phs.Cre.U = integer())
col <- c("Slc34a1_mRNA", "Scl34a3_mRNA", "Slc34a1_protein", "Slc34a3_protein", "Phs.Cre.U")

for (r in 1:length(col)){
  for (c in 1:length(col)){
      pval[r,c] <- cor.test(df6[,(4+r)], df6[,(4+c)], alternative = "two.sided", method = "pearson")$p.value
  }
}
rownames(pval) <- col
write.csv(pval, file = "./SpecificQ/Slc34a13_PhsCre_pval6.csv", row.names = FALSE, quote = FALSE)

# 12 mo.
df12 <- df[df$Age == 12,]
pval <- data.frame( Slc34a1_mRNA = integer(),
                    Scl34a3_mRNA = integer(),
                    Slc34a1_protein = integer(),
                    Slc34a3_protein = integer(),
                    Phs.Cre.U = integer())
col <- c("Slc34a1_mRNA", "Scl34a3_mRNA", "Slc34a1_protein", "Slc34a3_protein", "Phs.Cre.U")

for (r in 1:length(col)){
  for (c in 1:length(col)){
      pval[r,c] <- cor.test(df12[,(4+r)], df12[,(4+c)], alternative = "two.sided", method = "pearson")$p.value
  }
}
rownames(pval) <- col
write.csv(pval, file = "./SpecificQ/Slc34a13_PhsCre_pval12.csv", row.names = FALSE, quote = FALSE)

# 18 mo.
df18 <- df[df$Age == 18,]
pval <- data.frame( Slc34a1_mRNA = integer(),
                    Scl34a3_mRNA = integer(),
                    Slc34a1_protein = integer(),
                    Slc34a3_protein = integer(),
                    Phs.Cre.U = integer())
col <- c("Slc34a1_mRNA", "Scl34a3_mRNA", "Slc34a1_protein", "Slc34a3_protein", "Phs.Cre.U")

for (r in 1:length(col)){
  for (c in 1:length(col)){
      pval[r,c] <- cor.test(df18[,(4+r)], df18[,(4+c)], alternative = "two.sided", method = "pearson")$p.value
  }
}
rownames(pval) <- col
write.csv(pval, file = "./SpecificQ/Slc34a13_PhsCre_pval18.csv", row.names = FALSE, quote = FALSE)
