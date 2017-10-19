library(SNPtools)

snp.file <- "http://cgd.jax.org/tools/SNPtools/Build38/sanger.snps.NCBI38.txt.gz"
mgi.file <- "http://cgd.jax.org/tools/SNPtools/MGI/MGI.20130305.sorted.txt.gz"

available.strains <- get.strains(snp.file)
strains <- available.strains[c(4, 8, 2,  15, 16, 10, 17, 18)]

# Discontinued after 10/19/17 meeting
