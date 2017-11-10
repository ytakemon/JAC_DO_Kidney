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

df <- annot.samples
df$mrna <- expr.mrna[,m.id]
df$prot <- expr.protein[,p.id]

pdf(paste0("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Plot/", gene.name, "_mRNA-Protein_scatter_age.pdf"), height = 6, width = 6)
ggplot(df, aes(x = mrna, y = prot, colour = as.factor(Age))) +
  geom_point(alpah = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle(paste0(gene.name, " mRNA v. Protein"))
dev.off()
