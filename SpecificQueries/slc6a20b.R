library(ggplot2)
load("/hpcdata/cgd/QTL_mapping/kidney_combined/DO192_kidney.RData")
load("/hpcdata/cgd/shinyapps/kidney/shiny_annotation.RData")

gene.name <- "Slc6a20b"
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

m.id <- other.ids(gene.name, "mRNA")[[1]]
p.id <- other.ids(gene.name, "protein")[[3]]

df <- pheno
df$mrna <- expr.mrna[,m.id]
df$prot <- expr.protein[,p.id]

pdf("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Plot/Slc6a20b_scatter_age.pdf", height = 6, width = 6)
ggplot(df, aes(x = mrna, y = prot, colour = Age)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("Slc6a20b mRNA v. Protein")
dev.off()
