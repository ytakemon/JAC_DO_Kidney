# R/3.4.1
library(dplyr)
library(stringr)
library(ggpubr)
library(ggsci)
library(reshape2)
library(gridExtra)
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
Samples_akt$NZO_hom <- (str_sub(Samples_akt$all, 1, 1) == "E") & (str_sub(Samples_akt$all, 2, 2) == "")
Samples_akt$NZO_het <- (str_sub(Samples_akt$all, 1, 1) == "E") | (str_sub(Samples_akt$all, 2, 2) == "E")
Samples_akt_NZO_hom <- Samples_akt[Samples_akt$NZO_hom == TRUE, ]
Samples_akt_NZO_het <- Samples_akt[(Samples_akt$NZO_hom == FALSE & Samples_akt$NZO_het == TRUE), ]
Samples_akt_else <- Samples_akt[Samples_akt$NZO_het == FALSE, ]
# Subset pheno by group
Upheno$Sex <- Samples_akt$Sex
Pheno_NZO_hom <- Upheno[rownames(Upheno) %in% rownames(Samples_akt_NZO_hom),]
Pheno_NZO_het <- Upheno[rownames(Upheno) %in% rownames(Samples_akt_NZO_het),]
Pheno_else <- Upheno[rownames(Upheno) %in% rownames(Samples_akt_else),]
Pheno_NZO_hom$Allele <- "NZO/NZO"
Pheno_NZO_het$Allele <- "NZO/Other"
Pheno_else$Allele <- "Other/Other"
Pheno <- rbind(Pheno_NZO_hom, Pheno_NZO_het, Pheno_else)
# Normalize to Creatinine
Pheno$mg.cr.6 <- log(Pheno$mg.u.6 / Pheno$cr.u.6)
Pheno$mg.cr.12 <- log(Pheno$mg.u.12 / Pheno$cr.u.12)
Pheno$mg.cr.18 <- log(Pheno$mg.u.18 / Pheno$cr.u.18)
Pheno$ma.cr.6 <- log(Pheno$ma.u.6 / Pheno$cr.u.6)
Pheno$ma.cr.12 <- log(Pheno$ma.u.12 / Pheno$cr.u.12)
Pheno$ma.cr.18 <- log(Pheno$ma.u.18 / Pheno$cr.u.18)
Pheno$phs.cr.6 <- log(Pheno$phs.u.6 / Pheno$cr.u.6)
Pheno$phs.cr.12 <- log(Pheno$phs.u.12 / Pheno$cr.u.12)
Pheno$phs.cr.18 <- log(Pheno$phs.u.18 / Pheno$cr.u.18)

PhenoM <- Pheno[Pheno$Sex == "M",]
PhenoF <- Pheno[Pheno$Sex == "F",]

# Plots ----------------------------------------------------------------------
# T-test NZO vs non-NZO : Mg/Cr 6mo.
nNZO_hom <- nrow(PhenoM[complete.cases(PhenoM$mg.cr.6) & PhenoM$Allele == "NZO/NZO",])
nNZO_het <- nrow(PhenoM[complete.cases(PhenoM$mg.cr.6) & PhenoM$Allele == "NZO/Other",])
nOther <- nrow(PhenoM[complete.cases(PhenoM$mg.cr.6) & PhenoM$Allele == "Other/Other",])
nNZO_homF <- nrow(PhenoF[complete.cases(PhenoF$mg.cr.6) & PhenoF$Allele == "NZO/NZO",])
nNZO_hetF <- nrow(PhenoF[complete.cases(PhenoF$mg.cr.6) & PhenoF$Allele == "NZO/Other",])
nOtherF <- nrow(PhenoF[complete.cases(PhenoF$mg.cr.6) & PhenoF$Allele == "Other/Other",])
my_comparisons <- list( c(1, 2), c(2,3), c(3,1) )
Mg6 <- ggplot(Pheno, aes(x = Allele, y = mg.cr.6, colour = Allele))+
    geom_boxplot()+
    theme_bw()+
    labs( title = "log(Mg/Cr ratio) at 6 months",
          subtitle = paste0("Males: NZO/NZO = ",nNZO_hom, ", NZO/Other = ", nNZO_het,", Other = ", nOther, ")", "\n",
                          "Females: NZO/NZO = ",nNZO_homF, ", NZO/Other = ", nNZO_hetF,", Other = ", nOtherF, ")"),
          y = "log(Mg/Cr ratio) 6 months",
          x = "Founder Alleles")+
    scale_y_continuous(breaks = seq(-3, 3, 1))+
    facet_grid(. ~ Sex)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(label.y = 3, method = "anova")+
    guides(colour = FALSE)+
    scale_color_aaas()

# T-test NZO vs non-NZO : Mg/Cr 12mo.
nNZO_hom <- nrow(PhenoM[complete.cases(PhenoM$mg.cr.12) & PhenoM$Allele == "NZO/NZO",])
nNZO_het <- nrow(PhenoM[complete.cases(PhenoM$mg.cr.12) & PhenoM$Allele == "NZO/Other",])
nOther <- nrow(PhenoM[complete.cases(PhenoM$mg.cr.12) & PhenoM$Allele == "Other/Other",])
nNZO_homF <- nrow(PhenoF[complete.cases(PhenoF$mg.cr.12) & PhenoF$Allele == "NZO/NZO",])
nNZO_hetF <- nrow(PhenoF[complete.cases(PhenoF$mg.cr.12) & PhenoF$Allele == "NZO/Other",])
nOtherF <- nrow(PhenoF[complete.cases(PhenoF$mg.cr.12) & PhenoF$Allele == "Other/Other",])
my_comparisons <- list( c(1, 2), c(2,3), c(3,1) )
Mg12 <- ggplot(Pheno, aes(x = Allele, y = mg.cr.12, colour = Allele))+
    geom_boxplot()+
    theme_bw()+
    labs( title = "log(Mg/Cr ratio) at 12 months",
          subtitle = paste0("Males: NZO/NZO = ",nNZO_hom, ", NZO/Other = ", nNZO_het,", Other = ", nOther, ")", "\n",
                          "Females: NZO/NZO = ",nNZO_homF, ", NZO/Other = ", nNZO_hetF,", Other = ", nOtherF, ")"),
          y = "log(Mg/Cr ratio) 12 months",
          x = "Founder Alleles")+
    scale_y_continuous(breaks = seq(-8, 6, 1))+
    facet_grid(. ~ Sex)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(label.y = 6, method = "anova")+
    guides(colour = FALSE)+
    scale_color_aaas()

# T-test NZO vs non-NZO : Mg/Cr 18mo.
nNZO_hom <- nrow(PhenoM[complete.cases(PhenoM$mg.cr.18) & PhenoM$Allele == "NZO/NZO",])
nNZO_het <- nrow(PhenoM[complete.cases(PhenoM$mg.cr.18) & PhenoM$Allele == "NZO/Other",])
nOther <- nrow(PhenoM[complete.cases(PhenoM$mg.cr.18) & PhenoM$Allele == "Other/Other",])
nNZO_homF <- nrow(PhenoF[complete.cases(PhenoF$mg.cr.18) & PhenoF$Allele == "NZO/NZO",])
nNZO_hetF <- nrow(PhenoF[complete.cases(PhenoF$mg.cr.18) & PhenoF$Allele == "NZO/Other",])
nOtherF <- nrow(PhenoF[complete.cases(PhenoF$mg.cr.18) & PhenoF$Allele == "Other/Other",])
my_comparisons <- list( c(1, 2), c(2,3), c(3,1) )
Mg18 <- ggplot(Pheno, aes(x = Allele, y = mg.cr.18, colour = Allele))+
    geom_boxplot()+
    theme_bw()+
    labs( title = "log(Mg/Cr ratio) at 18 months",
          subtitle = paste0("Males: NZO/NZO = ",nNZO_hom, ", NZO/Other = ", nNZO_het,", Other = ", nOther, ")", "\n",
                          "Females: NZO/NZO = ",nNZO_homF, ", NZO/Other = ", nNZO_hetF,", Other = ", nOtherF, ")"),
          y = "log(Mg/Cr ratio) 18 months",
          x = "Founder Alleles")+
    scale_y_continuous(breaks = seq(-8, 6, 1))+
    facet_grid(. ~ Sex)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(label.y = 6, method = "anova")+
    guides(colour = FALSE)+
    scale_color_aaas()

# T-test NZO vs non-NZO : Alb/Cr 6mo. -----------------------------------------
nNZO_hom <- nrow(PhenoM[complete.cases(PhenoM$ma.cr.6) & PhenoM$Allele == "NZO/NZO",])
nNZO_het <- nrow(PhenoM[complete.cases(PhenoM$ma.cr.6) & PhenoM$Allele == "NZO/Other",])
nOther <- nrow(PhenoM[complete.cases(PhenoM$ma.cr.6) & PhenoM$Allele == "Other/Other",])
nNZO_homF <- nrow(PhenoF[complete.cases(PhenoF$ma.cr.6) & PhenoF$Allele == "NZO/NZO",])
nNZO_hetF <- nrow(PhenoF[complete.cases(PhenoF$ma.cr.6) & PhenoF$Allele == "NZO/Other",])
nOtherF <- nrow(PhenoF[complete.cases(PhenoF$ma.cr.6) & PhenoF$Allele == "Other/Other",])
my_comparisons <- list( c(2,3))
Alb6 <- ggplot(Pheno, aes(x = Allele, y = ma.cr.6, colour = Allele))+
    geom_boxplot()+
    theme_bw()+
    labs( title = "log(Alb/Cr ratio) at 6 months",
          subtitle = paste0("Males: NZO/NZO = ",nNZO_hom, ", NZO/Other = ", nNZO_het,", Other = ", nOther, ")", "\n",
                          "Females: NZO/NZO = ",nNZO_homF, ", NZO/Other = ", nNZO_hetF,", Other = ", nOtherF, ")"),
          y = "log(Alb/Cr ratio) 6 months",
          x = "Founder Alleles")+
    scale_y_continuous(breaks = seq(-8, 3, 1))+
    facet_grid(. ~ Sex)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(label.y = 3, method = "anova")+
    guides(colour = FALSE)+
    scale_color_aaas()

# T-test NZO vs non-NZO : Alb/Cr 12mo.
nNZO_hom <- nrow(PhenoM[complete.cases(PhenoM$ma.cr.12) & PhenoM$Allele == "NZO/NZO",])
nNZO_het <- nrow(PhenoM[complete.cases(PhenoM$ma.cr.12) & PhenoM$Allele == "NZO/Other",])
nOther <- nrow(PhenoM[complete.cases(PhenoM$ma.cr.12) & PhenoM$Allele == "Other/Other",])
nNZO_homF <- nrow(PhenoF[complete.cases(PhenoF$ma.cr.12) & PhenoF$Allele == "NZO/NZO",])
nNZO_hetF <- nrow(PhenoF[complete.cases(PhenoF$ma.cr.12) & PhenoF$Allele == "NZO/Other",])
nOtherF <- nrow(PhenoF[complete.cases(PhenoF$ma.cr.12) & PhenoF$Allele == "Other/Other",])
my_comparisons <- list( c(1, 2), c(2,3), c(3,1) )
Alb12 <- ggplot(Pheno, aes(x = Allele, y = ma.cr.12, colour = Allele))+
    geom_boxplot()+
    theme_bw()+
    labs( title = "log(Alb/Cr ratio) at 12 months",
          subtitle = paste0("Males: NZO/NZO = ",nNZO_hom, ", NZO/Other = ", nNZO_het,", Other = ", nOther, ")", "\n",
                          "Females: NZO/NZO = ",nNZO_homF, ", NZO/Other = ", nNZO_hetF,", Other = ", nOtherF, ")"),
          y = "log(Alb/Cr ratio) 12 months",
          x = "Founder Alleles")+
    scale_y_continuous(breaks = seq(-8, 6, 1))+
    facet_grid(. ~ Sex)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(label.y = 5, method = "anova")+
    guides(colour = FALSE)+
    scale_color_aaas()

# T-test NZO vs non-NZO : Alb/Cr 18mo.
nNZO_hom <- nrow(PhenoM[complete.cases(PhenoM$ma.cr.18) & PhenoM$Allele == "NZO/NZO",])
nNZO_het <- nrow(PhenoM[complete.cases(PhenoM$ma.cr.18) & PhenoM$Allele == "NZO/Other",])
nOther <- nrow(PhenoM[complete.cases(PhenoM$ma.cr.18) & PhenoM$Allele == "Other/Other",])
nNZO_homF <- nrow(PhenoF[complete.cases(PhenoF$ma.cr.18) & PhenoF$Allele == "NZO/NZO",])
nNZO_hetF <- nrow(PhenoF[complete.cases(PhenoF$ma.cr.18) & PhenoF$Allele == "NZO/Other",])
nOtherF <- nrow(PhenoF[complete.cases(PhenoF$ma.cr.18) & PhenoF$Allele == "Other/Other",])
my_comparisons <- list( c(2,3))
Alb18 <- ggplot(Pheno, aes(x = Allele, y = ma.cr.18, colour = Allele))+
    geom_boxplot()+
    theme_bw()+
    labs( title = "log(Alb/Cr ratio) at 18 months",
          subtitle = paste0("Males: NZO/NZO = ",nNZO_hom, ", NZO/Other = ", nNZO_het,", Other = ", nOther, ")", "\n",
                          "Females: NZO/NZO = ",nNZO_homF, ", NZO/Other = ", nNZO_hetF,", Other = ", nOtherF, ")"),
          y = "log(Alb/Cr ratio) 18 months",
          x = "Founder Alleles")+
    scale_y_continuous(breaks = seq(-8, 6, 1))+
    facet_grid(. ~ Sex)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(label.y = 3, method = "anova")+
    guides(colour = FALSE)+
    scale_color_aaas()

# T-test NZO vs non-NZO : Phs/Cr 6mo. -----------------------------------------
nNZO_hom <- nrow(PhenoM[complete.cases(PhenoM$phs.cr.6) & PhenoM$Allele == "NZO/NZO",])
nNZO_het <- nrow(PhenoM[complete.cases(PhenoM$phs.cr.6) & PhenoM$Allele == "NZO/Other",])
nOther <- nrow(PhenoM[complete.cases(PhenoM$phs.cr.6) & PhenoM$Allele == "Other/Other",])
nNZO_homF <- nrow(PhenoF[complete.cases(PhenoF$phs.cr.6) & PhenoF$Allele == "NZO/NZO",])
nNZO_hetF <- nrow(PhenoF[complete.cases(PhenoF$phs.cr.6) & PhenoF$Allele == "NZO/Other",])
nOtherF <- nrow(PhenoF[complete.cases(PhenoF$phs.cr.6) & PhenoF$Allele == "Other/Other",])
my_comparisons <- list( c(1, 2), c(2,3), c(3,1) )
Phs6 <- ggplot(Pheno, aes(x = Allele, y = phs.cr.6, colour = Allele))+
    geom_boxplot()+
    theme_bw()+
    labs( title = "log(Phs/Cr ratio) at 6 months",
          subtitle = paste0("Males: NZO/NZO = ",nNZO_hom, ", NZO/Other = ", nNZO_het,", Other = ", nOther, ")", "\n",
                          "Females: NZO/NZO = ",nNZO_homF, ", NZO/Other = ", nNZO_hetF,", Other = ", nOtherF, ")"),
          y = "log(Phs/Cr ratio) 6 months",
          x = "Founder Alleles")+
    scale_y_continuous(breaks = seq(-3, 6, 1))+
    facet_grid(. ~ Sex)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(label.y = 6, method = "anova")+
    guides(colour = FALSE)+
    scale_color_aaas()

# T-test NZO vs non-NZO : Phs/Cr 12mo.
nNZO_hom <- nrow(PhenoM[complete.cases(PhenoM$phs.cr.12) & PhenoM$Allele == "NZO/NZO",])
nNZO_het <- nrow(PhenoM[complete.cases(PhenoM$phs.cr.12) & PhenoM$Allele == "NZO/Other",])
nOther <- nrow(PhenoM[complete.cases(PhenoM$phs.cr.12) & PhenoM$Allele == "Other/Other",])
nNZO_homF <- nrow(PhenoF[complete.cases(PhenoF$phs.cr.12) & PhenoF$Allele == "NZO/NZO",])
nNZO_hetF <- nrow(PhenoF[complete.cases(PhenoF$phs.cr.12) & PhenoF$Allele == "NZO/Other",])
nOtherF <- nrow(PhenoF[complete.cases(PhenoF$phs.cr.12) & PhenoF$Allele == "Other/Other",])
my_comparisons <- list( c(1, 2), c(2,3), c(3,1) )
Phs12 <- ggplot(Pheno, aes(x = Allele, y = phs.cr.12, colour = Allele))+
    geom_boxplot()+
    theme_bw()+
    labs( title = "log(Phs/Cr ratio) at 12 months",
          subtitle = paste0("Males: NZO/NZO = ",nNZO_hom, ", NZO/Other = ", nNZO_het,", Other = ", nOther, ")", "\n",
                          "Females: NZO/NZO = ",nNZO_homF, ", NZO/Other = ", nNZO_hetF,", Other = ", nOtherF, ")"),
          y = "log(Phs/Cr ratio) 12 months",
          x = "Founder Alleles")+
    scale_y_continuous(breaks = seq(-8, 6, 1))+
    facet_grid(. ~ Sex)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(label.y = 6, method = "anova")+
    guides(colour = FALSE)+
    scale_color_aaas()

# T-test NZO vs non-NZO : Phs/Cr 18mo.
nNZO_hom <- nrow(PhenoM[complete.cases(PhenoM$phs.cr.18) & PhenoM$Allele == "NZO/NZO",])
nNZO_het <- nrow(PhenoM[complete.cases(PhenoM$phs.cr.18) & PhenoM$Allele == "NZO/Other",])
nOther <- nrow(PhenoM[complete.cases(PhenoM$phs.cr.18) & PhenoM$Allele == "Other/Other",])
nNZO_homF <- nrow(PhenoF[complete.cases(PhenoF$phs.cr.18) & PhenoF$Allele == "NZO/NZO",])
nNZO_hetF <- nrow(PhenoF[complete.cases(PhenoF$phs.cr.18) & PhenoF$Allele == "NZO/Other",])
nOtherF <- nrow(PhenoF[complete.cases(PhenoF$phs.cr.18) & PhenoF$Allele == "Other/Other",])
my_comparisons <- list( c(1, 2), c(2,3), c(3,1) )
Phs18 <- ggplot(Pheno, aes(x = Allele, y = phs.cr.18, colour = Allele))+
    geom_boxplot()+
    theme_bw()+
    labs( title = "log(Phs/Cr ratio) at 18 months",
          subtitle = paste0("Males: NZO/NZO = ",nNZO_hom, ", NZO/Other = ", nNZO_het,", Other = ", nOther, ")", "\n",
                          "Females: NZO/NZO = ",nNZO_homF, ", NZO/Other = ", nNZO_hetF,", Other = ", nOtherF, ")"),
          y = "log(Phs/Cr ratio) 18 months",
          x = "Founder Alleles")+
    scale_y_continuous(breaks = seq(-8, 6, 1))+
    facet_grid(. ~ Sex)+
    stat_compare_means(comparisons = my_comparisons)+
    stat_compare_means(label.y = 6, method = "anova")+
    guides(colour = FALSE)+
    scale_color_aaas()

pdf("./Plot/AktAllele_PhenoCompare.pdf", width = 20, height = 18)
grid.arrange(Mg6, Mg12, Mg18, Alb6, Alb12, Alb18, Phs6, Phs12, Phs18, ncol = 3,
             top = "Months", left = "Phenotype")
dev.off()

# Through time?
Pheno_time <- Pheno
Pheno_time <- Pheno_time[c(1,39:49)]
Pheno_time_Alb <- Pheno_time[,c(1,2,3,7:9)]
Pheno_time_Mg <- Pheno_time[,c(1,2,3,4:6)]
Pheno_time_Phs <- Pheno_time[,c(1,2,3,10:12)]
names(Pheno_time_Alb)[4:6] <- c("6mo", "12mo", "18mo")
names(Pheno_time_Mg)[4:6] <- c("6mo", "12mo", "18mo")
names(Pheno_time_Phs)[4:6] <- c("6mo", "12mo", "18mo")
Pheno_time_Alb  <- melt(Pheno_time_Alb)
Pheno_time_Mg  <- melt(Pheno_time_Mg)
Pheno_time_Phs  <- melt(Pheno_time_Phs)
names(Pheno_time_Alb)[4:5] <- c("TimePoint", "value")
names(Pheno_time_Mg)[4:5] <- c("TimePoint", "value")
names(Pheno_time_Phs)[4:5] <- c("TimePoint", "value")

Alb <- ggplot(Pheno_time_Alb, aes( x = TimePoint, y = value, colour = Allele))+
      geom_smooth(method = "lm", se=FALSE, aes(group = Allele, colour = Allele)) +
      geom_point(position = position_jitterdodge(), alpha = 0.5)+
      theme_bw()+
      facet_grid(. ~ Sex)+
      labs(title = "Log(Alb/Cr ratio) by Allele and Sex",
           y = "Log(Alb/Cr ratio)",
           x = "Time Point")+
      guides( colour = FALSE)+
      scale_color_aaas()

Mg <- ggplot(Pheno_time_Mg, aes( x = TimePoint, y = value, colour = Allele))+
      geom_smooth(method = "lm", se=FALSE, aes(group = Allele, colour = Allele)) +
      geom_point(position = position_jitterdodge(), alpha = 0.5)+
      theme_bw()+
      labs(title = "Log(Mg/Cr ratio) by Allele and Sex",
           y = "Log(Mg/Cr ratio)",
           x = "Time Point")+
      facet_grid(. ~ Sex)+
      guides( colour = FALSE)+
      scale_color_aaas()

Phs <- ggplot(Pheno_time_Phs, aes( x = TimePoint, y = value, colour = Allele))+
      geom_smooth(method = "lm", se=FALSE, aes(group = Allele, colour = Allele)) +
      geom_point(position = position_jitterdodge(), alpha = 0.5)+
      theme_bw()+
      labs(title = "Log(Phs/Cr ratio) by Allele and Sex",
           y = "Log(Phs/Cr ratio)",
           x = "Time Point")+
      facet_grid(. ~ Sex)+
      scale_color_aaas()

pdf("./Plot/AktAllele_PhenoCompare_time.pdf", width = 20, height = 6)
grid.arrange(Mg, Alb, Phs, ncol = 3)
dev.off()
