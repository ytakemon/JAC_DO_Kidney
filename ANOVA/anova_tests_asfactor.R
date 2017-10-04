# ANOVA
#
# Calculates averages and fit several ANOVA models.
#
# References:
# https://stats.stackexchange.com/questions/76815/multiple-regression-or-partial-correlation-coefficient-and-relations-between-th
# https://stats.stackexchange.com/questions/32464/how-does-the-correlation-coefficient-differ-from-regression-slope
#
# Usage:
#
# Rscript anova_tests.R input.RData output.csv

library(broom)
# library(ppcor) # library for partial correlation coefs

# Options -----------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# for kidney
# args <- list("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/RNAseq_data/DO188b_kidney_noprobs.RData", "/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Anova_output/kidney_anova_factor_table.csv")
#

# two arguments expected (input and output files)
stopifnot(length(args)==2)
input.file = args[[1]]
stopifnot(file.exists(input.file))
output.file = args[[2]]
load(input.file)

# ANOVA FACTORS (p-values, normalized coefs, and slopes) -----------------------------------

## test for dependence between Age/Sex as factors and x (x is mRNA/Prot expression)
anova_tests_2 <- function(x) {
  # Find adjusted p.values using TukeyHSD
  fit_aov <- aov(x ~ Age + Sex + Generation, data=annot.samples)
  HSD <- TukeyHSD(fit_aov, ordered = TRUE)$Age
  # Determine direction for HSDx
  if (any(rownames(HSD) == "6-12") | any(rownames(HSD) == "12-6")){
    if(any(rownames(HSD) == "6-12")){
      age6_12 <- "6-12"
    } else if (any(rownames(HSD) == "12-6")){
      age6_12 <- "12-6"
    }
  }
  if (any(rownames(HSD) == "12-18") | any(rownames(HSD) == "18-12")){
    if(any(rownames(HSD) == "12-18")){
      age12_18 <- "12-18"
    } else if (any(rownames(HSD) == "18-12")){
      age12_18 <- "18-12"
    }
  }
  if (any(rownames(HSD) == "6-18") | any(rownames(HSD) == "18-6")){
    if(any(rownames(HSD) == "6-18")){
      age6_18 <- "6-18"
    } else if (any(rownames(HSD) == "18-6")){
      age6_18 <- "18-6"
    }
  }
  pvalues <- c(HSD[age6_12, "p adj"],
               HSD[age12_18, "p adj"],
               HSD[age6_18, "p adj"])
  # Find coefficients via individual lm per comparison group
  # 6-12
  sub <- which(annot.samples$Age %in% c(6, 12))
  y <- x[sub]
  fit_lm <- lm(y ~ Age + Sex + Generation, data=annot.samples[annot.samples$Age %in% c(6,12),])
  slope6_12 <- fit_lm$coefficients[["Age12"]]
 # 12-18
  sub <- which(annot.samples$Age %in% c(12, 18))
  y <- x[sub]
  fit_lm <- lm(y ~ Age + Sex + Generation, data=annot.samples[annot.samples$Age %in% c(12,18),])
  slope12_18 <- fit_lm$coefficients[["Age18"]]
 # 6-18
  sub <- which(annot.samples$Age %in% c(6, 18))
  y <- x[sub]
  fit_lm <- lm(y ~ Age + Sex + Generation, data=annot.samples[annot.samples$Age %in% c(6,18),])
  slope6_18 <- fit_lm$coefficients[["Age18"]]
  slopes <- c(slope6_12, slope12_18, slope6_18)

  return(c(pvalues,slopes))
}

## test for interaction between Age and Sex
anova_tests_int <- function(x,y) {
  # full model with interaction
  fitx_aov <- aov(x ~ Age * Sex + Generation, data=annot.samples)
  fity_aov <- aov(y ~ Age * Sex + Generation, data=annot.samples)
  HSDx <- TukeyHSD(fitx_aov, which = "Age", ordered = TRUE)$Age
  HSDy <- TukeyHSD(fity_aov, which = "Age", ordered = TRUE)$Age # weirdly does not like to flip names, but values have changed
  # Determine direction for HSDx
  if (any(rownames(HSDx) == "6-12") | any(rownames(HSDx) == "12-6")){
    if(any(rownames(HSDx) == "6-12")){
      agex6_12 <- "6-12"
    } else if (any(rownames(HSDx) == "12-6")){
      agex6_12 <- "12-6"
    }
  }
  if (any(rownames(HSDx) == "12-18") | any(rownames(HSDx) == "18-12")){
    if(any(rownames(HSDx) == "12-18")){
      agex12_18 <- "12-18"
    } else if (any(rownames(HSDx) == "18-12")){
      agex12_18 <- "18-12"
    }
  }
  if (any(rownames(HSDx) == "6-18") | any(rownames(HSDx) == "18-6")){
    if(any(rownames(HSDx) == "6-18")){
      agex6_18 <- "6-18"
    } else if (any(rownames(HSDx) == "18-6")){
      agex6_18 <- "18-6"
    }
  }
  # Determine direction for HSDy
  if (any(rownames(HSDy) == "6-12") | any(rownames(HSDy) == "12-6")){
    if(any(rownames(HSDy) == "6-12")){
      agey6_12 <- "6-12"
    } else if (any(rownames(HSDy) == "12-6")){
      agey6_12 <- "12-6"
    }
  }
  if (any(rownames(HSDy) == "12-18") | any(rownames(HSDy) == "18-12")){
    if(any(rownames(HSDy) == "12-18")){
      agey12_18 <- "12-18"
    } else if (any(rownames(HSDy) == "18-12")){
      agey12_18 <- "18-12"
    }
  }
  if (any(rownames(HSDy) == "6-18") | any(rownames(HSDy) == "18-6")){
    if(any(rownames(HSDy) == "6-18")){
      agey6_18 <- "6-18"
    } else if (any(rownames(HSDy) == "18-6")){
      agey6_18 <- "18-6"
    }
  }

  # Determine p-value
  pvalues <- c(HSDx[agex6_12, "p adj"],
               HSDx[agex12_18, "p adj"],
               HSDx[agex6_18, "p adj"],
               HSDy[agey6_12, "p adj"],
               HSDy[agey12_18, "p adj"],
               HSDy[agey6_18, "p adj"])

  # Find coefficients via individual lm per comparison group
  # x = mRNA
  # 6-12
  sub <- which(annot.samples$Age %in% c(6, 12))
  x_sub <- x[sub]
  fit_lm <- lm(x_sub ~ Age * Sex + Generation, data=annot.samples[annot.samples$Age %in% c(6,12),])
  slopex6_12 <- fit_lm$coefficients[["Age12"]]
 # 12-18
  sub <- which(annot.samples$Age %in% c(12, 18))
  x_sub <- x[sub]
  fit_lm <- lm(x_sub ~ Age * Sex + Generation, data=annot.samples[annot.samples$Age %in% c(12,18),])
  slopex12_18 <- fit_lm$coefficients[["Age18"]]
 # 6-18
  sub <- which(annot.samples$Age %in% c(6, 18))
  x_sub <- x[sub]
  fit_lm <- lm(x_sub ~ Age * Sex + Generation, data=annot.samples[annot.samples$Age %in% c(6,18),])
  slopex6_18 <- fit_lm$coefficients[["Age18"]]
  slopex <- c(slopex6_12, slopex12_18, slopex6_18)

  # y = Protein
  # 6-12
  sub <- which(annot.samples$Age %in% c(6, 12))
  y_sub <- y[sub]
  fit_lm <- lm(y_sub ~ Age * Sex + Generation, data=annot.samples[annot.samples$Age %in% c(6,12),])
  slopey6_12 <- fit_lm$coefficients[["Age12"]]
 # 12-18
  sub <- which(annot.samples$Age %in% c(12, 18))
  y_sub <- y[sub]
  fit_lm <- lm(y_sub ~ Age * Sex + Generation, data=annot.samples[annot.samples$Age %in% c(12,18),])
  slopey12_18 <- fit_lm$coefficients[["Age18"]]
 # 6-18
  sub <- which(annot.samples$Age %in% c(6, 18))
  y_sub <- y[sub]
  fit_lm <- lm(y_sub ~ Age * Sex + Generation, data=annot.samples[annot.samples$Age %in% c(6,18),])
  slopey6_18 <- fit_lm$coefficients[["Age18"]]
  slopey <- c(slopey6_12, slopey12_18, slopey6_18)

  return(c(pvalues, slopex, slopey))
}

cols <- c("p.mRNA_Int_fact_6_12", "p.mRNA_Int_fact_12_18", "p.mRNA_Int_fact_6_18",
          "p.Prot_Int_fact_6_12", "p.Prot_Int_fact_12_18", "p.Prot_Int_fact_6_18",
          "m.mRNA_Int_fact_6_12", "m.mRNA_Int_fact_12_18", "m.mRNA_Int_fact_6_18",
          "m.Prot_Int_fact_6_12", "m.Prot_Int_fact_12_18", "m.Prot_Int_fact_6_18",
          "p.mRNA_Age.Sex_fact_6_12","p.mRNA_Age.Sex_fact_12_18","p.mRNA_Age.Sex_fact_6_18",
          "m.mRNA_Age.Sex_fact_6_12","m.mRNA_Age.Sex_fact_12_18","m.mRNA_Age.Sex_fact_6_18",
          "p.Prot_Age.Sex_fact_6_12","p.Prot_Age.Sex_fact_12_18","p.Prot_Age.Sex_fact_6_18",
          "m.Prot_Age.Sex_fact_6_12","m.Prot_Age.Sex_fact_12_18","m.Prot_Age.Sex_fact_6_18")

result.table <- matrix(NA, N[["pairs"]], length(cols))
colnames(result.table) <- cols
annot.samples$Age <- factor(annot.samples$Age)
annot.samples$Sex <- factor(annot.samples$Sex)

## normalize verything to mean=0, sd=1
## for easy calculation of partial correlation coefficients
for (i in 1:N[["pairs"]]) {
  expr.mrna[,i] <- scale(expr.mrna[,i])
  expr.protein[,i] <- scale(expr.protein[,i])
}

# interaction tests as a factor
print("Testing for interaction between Age and Sex as factors...")
for (i in 1:N[["pairs"]]) {
  if (i %% 100 == 0) print(i)
  cols <- c("p.mRNA_Int_fact_6_12", "p.mRNA_Int_fact_12_18", "p.mRNA_Int_fact_6_18",
            "p.Prot_Int_fact_6_12", "p.Prot_Int_fact_12_18", "p.Prot_Int_fact_6_18",
            "m.mRNA_Int_fact_6_12", "m.mRNA_Int_fact_12_18", "m.mRNA_Int_fact_6_18",
            "m.Prot_Int_fact_6_12", "m.Prot_Int_fact_12_18", "m.Prot_Int_fact_6_18")
  result.table[i, which(colnames(result.table) %in% cols)] <- anova_tests_int(expr.mrna[,i], expr.protein[,i])
}

# ANOVA for mRNA as a factor
print("Testing for dependence between mRNA and Age/Sex as factors...")
for (i in 1:N[["pairs"]]) {
  if (i %% 100 == 0) print(i)
  result.table[i, which(colnames(result.table) %in% c("p.mRNA_Age.Sex_fact_6_12","p.mRNA_Age.Sex_fact_12_18","p.mRNA_Age.Sex_fact_6_18",
                                                      "m.mRNA_Age.Sex_fact_6_12","m.mRNA_Age.Sex_fact_12_18","m.mRNA_Age.Sex_fact_6_18"))] <- anova_tests_2(expr.mrna[,i])
}

# ANOVA for proteins as a factor
print("Testing for dependence between mRNA and Age/Sex as factors...")
for (i in 1:N[["pairs"]]) {
  if (i %% 100 == 0) print(i)
  result.table[i, which(colnames(result.table) %in% c("p.Prot_Age.Sex_fact_6_12","p.Prot_Age.Sex_fact_12_18","p.Prot_Age.Sex_fact_6_18",
                                                      "m.Prot_Age.Sex_fact_6_12","m.Prot_Age.Sex_fact_12_18","m.Prot_Age.Sex_fact_6_18"))] <- anova_tests_2(expr.protein[,i])
}

# reorder columns - p-values first, corelation coefs second
pcols <- grep("^p.", colnames(result.table))
mcols <- grep("^m.", colnames(result.table))
result.table <- result.table[,c(pcols,mcols)]

# Output ------------------------------------------------------------------

annot.cols <- c("id","gene_id","symbol","chr","start","end","strand","biotype")

output <- cbind(annot.protein[1:N[["pairs"]], annot.cols],
                result.table)

write.csv(output, output.file, row.names=FALSE)
