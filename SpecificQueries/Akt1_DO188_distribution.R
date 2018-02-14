# R/3.4.1
library(dplyr)
library(stringr)
library(ggpubr)
library(ggsci)
library(reshape2)
library(gridExtra)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")

# Annotate Akt1 founder allele
# Get snps markers
# Akt1 pos : 12: 112,653,821-112,674,884
Akt_snps <- snps[snps$chr == "12",]
Akt_snps <- Akt_snps[(Akt_snps$bp > 112.6 & Akt_snps$bp < 112.7),]
# Get probs for all markers
probs <- genoprobs[,, colnames(genoprobs[1,,]) %in% Akt_snps$marker]
#> dim(probs)
#[1] 188    8    3

# Get founder allele that has probs > 0.49
# First check:
identical(rownames(annot.samples), rownames(probs))
#> dim(annot.samples)
#[1] 1045   17
Samples_akt <- annot.samples
Samples_akt$snp1 <- NA
Samples_akt$snp2 <- NA
Samples_akt$snp3 <- NA

# Annotate 3 SNPs
for ( i in 1:length(rownames(Samples_akt))){
  for ( n in 1:3){
    temp <- probs[i,,n]
    temp <- names(temp[temp > 0.30])

    if(length(temp) == 2){
      temp <- paste0(temp[1], temp[2])
    } else if (length(temp == 1)){
      temp <- temp
    } else if (length(temp > 2)){
      print("Error: more than 2 founders!")
    }
    Samples_akt[i ,(10+n)] <- temp
  }
}

# Get most likely genotype
Samples_akt$all <- NA
for (i in 1:length(rownames(Samples_akt))){
  Samples_akt$all[i] <- max(t(Samples_akt[i,11:13]))
}
# NZO allele TRUE/FALSE
Samples_akt$NZO_hom <- (str_sub(Samples_akt$all, 1, 1) == "E") & (str_sub(Samples_akt$all, 2, 2) == "")
Samples_akt$NZO_het <- (str_sub(Samples_akt$all, 1, 1) == "E") | (str_sub(Samples_akt$all, 2, 2) == "E")

#Assign geno: "NZO/NZO", "NZO/Other", "Other/Other"
Samples_akt$geno <- NA
Samples_akt[Samples_akt$NZO_hom == TRUE, ]$geno <- "NZO/NZO"
Samples_akt[(Samples_akt$NZO_hom == FALSE & Samples_akt$NZO_het == TRUE), ]$geno <- "NZO/Other"
Samples_akt_else <- Samples_akt[Samples_akt$NZO_het == FALSE, ]$geno <- "Other/Other"

# Add Akt1 mRNA expression
id <- annot.mrna[annot.mrna$symbol == "Akt1",]$id
pid <- annot.protein[annot.protein$symbol == "Akt1",]$id
Samples_akt$Akt1 <- expr.mrna[,id]
Samples_akt$Akt1_raw <- raw.mrna[,id]
Samples_akt$Akt1_prot <- expr.protein[,pid]
Samples_akt$Akt1_prot_raw <- raw.protein[,pid]

# Age needs to be a factor
Samples_akt$Age <- as.factor(as.character(Samples_akt$Age))
Samples_akt$Age <- factor(Samples_akt$Age, levels(Samples_akt$Age)[c(3,1,2)])

# Plot mRNA
pdf("./Plot/Akt1Allele_Akt1exp_by_Age.pdf", width = 12, height = 7)
ggplot(Samples_akt, aes(x = Age, y = Akt1_raw, colour = geno)) +
  geom_smooth(method = "lm", se = FALSE, aes(group = geno, colour = geno)) +
  theme_bw() +
  labs( title = "Akt1 mRNA expression by genotype",
        subtitle = paste0("Female: 93, (NZO/NZO=0 , NZO/Other=9 , Other/Other=84 ) \nMale: 95, (NZO/NZO=3 , NZO/Other=12 , Other/Other=80 )"),
        y = "Akt1 mRNA expression",
        x = "Age") +
  facet_grid(. ~ Sex) +
  scale_colour_aaas()
dev.off()

# Plot protein
pdf("./Plot/Akt1Allele_Akt1prot_by_Age.pdf", width = 12, height = 7)
ggplot(Samples_akt, aes(x = Age, y = Akt1_prot_raw, colour = geno)) +
  geom_smooth(method = "lm", se = FALSE, aes(group = geno, colour = geno)) +
  theme_bw() +
  labs( title = "Akt1 protein expression by genotype",
        subtitle = paste0("Female: 93, (NZO/NZO=0 , NZO/Other=9 , Other/Other=84 ) \nMale: 95, (NZO/NZO=3 , NZO/Other=12 , Other/Other=80 )"),
        y = "Akt1 protein expression",
        x = "Age") +
  facet_grid(. ~ Sex) +
  scale_colour_aaas()
dev.off()
