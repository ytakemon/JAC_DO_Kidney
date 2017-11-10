library(ggplot2)
library(dplyr)
library(gridExtra)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")

# Identify gene name
gene_name <- "Pdgfb"
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
gene <- other.ids(gene_name, "mRNA")

# new df
df <- annot.samples[,1:4]
df$mRNA <- expr.mrna[, gene$id]
df$Protein <- expr.protein[, gene$protein_id]
df$Age <- as.factor(df$Age)
df$Sex <- as.factor(df$Sex)


mRNA <- ggplot(df, aes(x = Age, y = mRNA, fill = Sex))+
              geom_smooth(method = "lm", se=FALSE, aes(group = Sex, colour = Sex, fill = Sex)) +
              geom_point(aes(colour = Sex, fill = Sex, alpha = 0.5), size = 2, shape = 21, position = position_jitterdodge())+
              labs( title = paste0(gene, " ", gene_name, " mRNA")) +
              guides(alpha = FALSE) +
              theme(panel.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(color = "black"))

protein <- ggplot(df, aes(x = Age, y = Protein, fill = Sex))+
              geom_smooth(method = "lm", se=FALSE, aes(group = Sex, colour = Sex, fill = Sex)) +
              geom_point(aes(colour = Sex, fill = Sex, alpha = 0.5), size = 2, shape = 21, position = position_jitterdodge())+
              labs(title = paste0(gene, " ", gene_name, " Protein")) +
              guides(alpha = FALSE) +
              theme(panel.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(color = "black"))

pdf(paste0("./Plot/", gene_name,"_Age.pdf"), width = 4, height = 3)
grid.arrange(mRNA, protein, ncol = 2)
dev.off()
