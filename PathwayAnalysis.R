library(ReactomePA)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
file <- "kidney_anova_slope_output.csv"
data <- read.csv(paste0("./Anova_output/",file), header = T)

# set up dan's function
rankZ = function(x) {
  x = rank(x, na.last = "keep", ties.method = "average") / (sum(!is.na(x)) + 1)
  return(qnorm(x))
}

# grab significant RNA
sig_mRNA <- data %>%
  filter(p.mRNA_Age.Sex < 0.05) %>%
  select(id, gene_id, symbol, chr, start, end, strand, biotype, p.mRNA_Age.Sex) %>%
  mutate(rankZ = rankZ(p.mRNA_Age.Sex),
         ENTREZID = NA)

# translate ensembl id to EntrezID
for(i in 1:nrow(sig_mRNA)){
  print(i)
  tryCatch({
    sig_mRNA$ENTREZID[i] <- bitr(as.character(sig_mRNA$gene_id[i]), fromType = "ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
  }, error = function(e){
    sig_mRNA$ENTREZID[i] <- NA 
  })
}

# no ENTREZ ID found for this many genes:
nrow(sig_mRNA[is.na(sig_mRNA$ENTREZID),])
head(sig_mRNA[is.na(sig_mRNA$ENTREZID),])



keytypes(org.Mm.eg.db)
sig_mRNA[is.na(sig_mRNA$ENTREZID),]

# source
https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
