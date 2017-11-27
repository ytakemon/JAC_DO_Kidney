library(dplyr)
library(stringr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO1045_kidney.Rdata")

# Number of samples
#> dim(Upheno)
#[1] 1045   38
#> dim(genoprobs)
#[1]  1045     8 68268
#> dim(samples)
#[1] 1045   17

# Annotate Akt1 founder allele

# Get snps markers
# Akt1 pos : 12: 112,653,821-112,674,884
Akt_snps <- MM_snps[MM_snps$chr == "12",]
Akt_snps <- Akt_snps[(Akt_snps$pos > 112.6 & Akt_snps$pos < 112.7),]
# Get probs for all markers
probs <- genoprobs[,, colnames(genoprobs[1,,]) %in% Akt_snps$marker]
#> dim(probs)
#[1] 1045    8    9

# Get founder allele that has probs > 0.49
# First check:
identical(rownames(samples), rownames(probs))
#> dim(samples)
#[1] 1045   17
Samples_akt <- samples
Samples_akt$snp1 <- NA
Samples_akt$snp2 <- NA
Samples_akt$snp3 <- NA
Samples_akt$snp4 <- NA
Samples_akt$snp5 <- NA

# Annotate 9 SNPs
for ( i in 1:length(rownames(Samples_akt))){
  for ( n in 1:5){
    temp <- probs[i,,n]
    temp <- names(temp[temp > 0.30])

    if(length(temp) == 2){
      temp <- paste0(temp[1], temp[2])
    } else if (length(temp == 1)){
      temp <- temp
    } else if (length(temp > 2)){
      print("Error: more than 2 founders!")
    }
    Samples_akt[i ,(17+n)] <- temp
  }
}

# Get most likely genotype
Samples_akt$all <- NA
for (i in 1:length(rownames(Samples_akt))){
  Samples_akt$all[i] <- max(t(Samples_akt[i,18:22]))
}
# NZO allele TRUE/FALSE
Samples_akt$NZO <- (str_sub(Samples_akt$all, 1, 1) == "E") | (str_sub(Samples_akt$all, 2, 2) == "E")
Samples_akt_NZO <- Samples_akt[Samples_akt$NZO == TRUE, ]
Samples_akt_else <- Samples_akt[Samples_akt$NZO == FALSE, ]
# Subset pheno by group
Pheno_NZO <- Upheno[rownames(Upheno) %in% rownames(Samples_akt_NZO),]
Pheno_else <- Upheno[rownames(Upheno) %in% rownames(Samples_akt_else),]
#> dim(Pheno_NZO)
#[1] 136  38
#> dim(Pheno_else)
#[1] 909  38

# T-test NZO vs non-NZO : Mg/Cr 6mo.
# T-test NZO vs non-NZO : Mg/Cr 12mo.
# T-test NZO vs non-NZO : Mg/Cr 18mo.

# T-test NZO vs non-NZO : Alb/Cr 6mo.
# T-test NZO vs non-NZO : Alb/Cr 12mo.
# T-test NZO vs non-NZO : Alb/Cr 18mo.

# T-test NZO vs non-NZO : Phs/Cr 6mo.
# T-test NZO vs non-NZO : Phs/Cr 12mo.
# T-test NZO vs non-NZO : Phs/Cr 18mo.
