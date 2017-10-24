library(ggplot2)
library(dplyr)
library(gridExtra)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO188b_kidney.RData")
pheno <- read.delim("./Phenotype/phenotypes/formatted/JAC_CS_urine_chem_v1.txt", sep = "\t")

# Subset pheno to match count for animal id
pheno$duplicated <- (duplicated(pheno$mouse.id) | duplicated(pheno$mouse.id, fromLast = TRUE))
pheno <- pheno[pheno$duplicated == FALSE,]
pheno <- pheno[pheno$mouse.id %in% annot.samples$Mouse.ID,] # only 176 animals with phenotypes
pheno <- arrange(pheno, mouse.id)
pheno$log.ma.cr.u <- log(pheno$ma.cr.u)
for( i in 1:length(pheno$log.ma.cr.u)){
  if(is.na(pheno$log.ma.cr.u[i])){
    pheno$log.ma.cr.u[i] <- NA
  } else if(pheno$log.ma.cr.u[i] == -Inf){
    pheno$log.ma.cr.u[i] <- 0
  } else if (pheno$log.ma.cr.u[i] < 0){
    pheno$log.ma.cr.u[i] <- NA
  } else {
    pheno$log.ma.cr.u[i] <- pheno$log.ma.cr.u[i]
  }
}

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
colnames(df)[5:9] <- c("Slc34a1_mRNA", "Slc34a3_mRNA", "Slc34a1_protein", "Slc34a3_protein", "Phs.Cre.U")
df$Age <- as.factor(df$Age)
df$Sex <- as.factor(df$Sex)



mRNA1 <- ggplot(df, aes(x = Age, y = Slc34a1_mRNA, fill = Sex))+
              geom_smooth(method = "lm", se=FALSE, aes(group = Sex, colour = Sex, fill = Sex)) +
              geom_point(aes(fill = Sex), size = 2, shape = 21, position = position_jitterdodge())+
              ggtitle("Slc34a1 mRNA") +
              theme(panel.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(color = "black"))

protein1 <- ggplot(df, aes(x = Age, y = Slc34a1_protein, fill = Sex))+
              geom_smooth(method = "lm", se=FALSE, aes(group = Sex, colour = Sex, fill = Sex)) +
              geom_point(aes(fill = Sex), size = 2, shape = 21, position = position_jitterdodge())+
              ggtitle("Slc34a1 protein") +
              theme(panel.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(color = "black"))

mRNA3 <- ggplot(df, aes(x = Age, y = Slc34a3_mRNA, fill = Sex))+
              geom_smooth(method = "lm", se=FALSE, aes(group = Sex, colour = Sex, fill = Sex)) +
              geom_point(aes(fill = Sex), size = 2, shape = 21, position = position_jitterdodge())+
              ggtitle("Slc34a3 mRNA") +
              theme(panel.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(color = "black"))

protein3 <- ggplot(df, aes(x = Age, y = Slc34a3_protein, fill = Sex))+
              geom_smooth(method = "lm", se=FALSE, aes(group = Sex, colour = Sex, fill = Sex)) +
              geom_point(aes(fill = Sex), size = 2, shape = 21, position = position_jitterdodge())+
              ggtitle("Slc34a3 protein") +
              theme(panel.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(color = "black"))

pdf("./Plot/Slc34a13_Age.pdf", width = 7, height = 6)
grid.arrange(mRNA1, protein1, mRNA3, protein3, ncol = 2)
dev.off()
