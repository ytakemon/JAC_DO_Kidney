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
# args <- list("DO188b_kidney_noprobs.RData", "mrna.kidney_anova_table.csv", "protein.kidney_anova_table.csv")

# two arguments expected (input and output files)
stopifnot(length(args)==3)

input.file = list.files(path = "/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/RNAseq_data", pattern = paste0("^",args[[1]]),recursive = TRUE,full.names = TRUE)
stopifnot(file.exists(input.file))
m.output.file = args[[2]]
p.output.file = args[[3]]
load(input.file)

# ANOVA (p-values, normalized coefs, and slopes) -----------------------------------
## test for dependence between Age/Sex and x (x is mRNA/Prot expression)
anova_tests_2 <- function(x) {
  # full model without interaction
  tmp.full <- lm(x ~ Age + Sex + Generation, data=annot.samples)
  lm.full <- tidy(tmp.full)

  pvalues <- c(subset(lm.full, term=="Age")$p.value,
               subset(lm.full, term=="Sex")$p.value)

  # partial correlation coefficient equals normalized beta * CONST
  # see http://stats.stackexchange.com/questions/76815/multiple-regression-or-partial-correlation-coefficient-and-relations-between-th

  pres <- !is.na(x) # must be calculated on the same observations
  tmp.age  <- lm(Age ~ Sex + Generation, data=annot.samples[pres,])
  tmp.sex  <- lm(Sex ~ Age + Generation, data=annot.samples[pres,])
  sigma.full <- sd(tmp.full$resid, na.rm = TRUE)
  sigma.age <- sd(tmp.age$resid, na.rm = TRUE)
  sigma.sex <- sd(tmp.sex$resid, na.rm = TRUE)

  coefs <- c(subset(lm.full, term=="Age")$estimate * sigma.age,
             subset(lm.full, term=="Sex")$estimate * sigma.sex) / sigma.full
  slopes <- c(subset(lm.full, term=="Age")$estimate,
             subset(lm.full, term=="Sex")$estimate)

  return(c(pvalues,coefs,slopes))
}

# Test for sex:age interaction
anova_tests_int <- function(x) {
  # full model with interaction
  lm.x <- tidy(lm(x ~ Age * Sex + Generation, data=annot.samples))
  pvalues <- subset(lm.x, term=="Age:Sex")$p.value

  return(pvalues)
}


# Remove duplicates in mrna objects, duplicates were previously created becase
# there were multiple annotated protein that corresponded to a transcript.
# For this analysis, we are not combining rna to protein, thus duplicates need
# to be removed
annot.mrna <- annot.mrna[which(annot.mrna$duplicated == FALSE),]
expr.mrna <- expr.mrna[,!duplicated(colnames(expr.mrna))]

mcols <- c("p.mRNA_Age.Sex", "p.mRNA_Sex.Age",
           "r.mRNA_Age.Sex", "r.mRNA_Sex.Age",
           "m.mRNA_Age.Sex", "m.mRNA_Sex.Age",
           "p.mRNA_AgeSexInt")

pcols <- c("p.Prot_Age.Sex", "p.Prot_Sex.Age",
           "r.Prot_Age.Sex", "r.Prot_Sex.Age",
           "m.Prot_Age.Sex", "m.Prot_Sex.Age",
           "p.Prot_AgeSexInt")


m.result.table <- matrix(NA, dim(expr.mrna)[2], length(mcols))
p.result.table <- matrix(NA, dim(expr.protein)[2], length(pcols))
colnames(m.result.table) <- mcols
colnames(p.result.table) <- pcols

## normalize verything to mean=0, sd=1
## for easy calculation of partial correlation coefficients
for (i in 1:dim(expr.mrna)[2]) {
  expr.mrna[,i] <- scale(expr.mrna[,i])
}
for (i in 1:dim(expr.protein)[2]) {
  expr.protein[,i] <- scale(expr.protein[,i])
}
annot.samples$Age <- scale(annot.samples$Age)
annot.samples$Sex <- scale(as.numeric(factor(annot.samples$Sex)))

# ANOVA for mRNA
print("Testing for dependence between mRNA and Age/Sex...")
for (i in 1:dim(expr.mrna)[2]) {
  if (i %% 100 == 0) print(i)
  m.result.table[i, which(colnames(m.result.table) %in% c("p.mRNA_Age.Sex","p.mRNA_Sex.Age",
    "r.mRNA_Age.Sex","r.mRNA_Sex.Age",
    "m.mRNA_Age.Sex","m.mRNA_Sex.Age"))] <- anova_tests_2(expr.mrna[,i])
}

print("Testing for interaction between mRNA and Age:Sex interaction...")
for (i in 1:dim(expr.mrna)[2]) {
  if (i %% 100 == 0) print(i)
  m.result.table[i, "p.mRNA_AgeSexInt"] <- anova_tests_int(expr.mrna[,i])
}

# ANOVA for proteins
print("Testing for dependence between protein and Age/Sex...")
for (i in 1:dim(expr.protein)[2]) {
  if (i %% 100 == 0) print(i)
  p.result.table[i, which(colnames(p.result.table) %in% c("p.Prot_Age.Sex","p.Prot_Sex.Age",
    "r.Prot_Age.Sex","r.Prot_Sex.Age",
    "m.Prot_Age.Sex","m.Prot_Sex.Age"))] <- anova_tests_2(expr.protein[,i])
}

print("Testing for interaction between protein and Age:Sex...")
for (i in 1:dim(expr.protein)[2]) {
  if (i %% 100 == 0) print(i)
  p.result.table[i, "p.Prot_AgeSexInt"] <- anova_tests_int(expr.protein[,i])
}


# reorder columns - p-values first, corelation coefs second

pcols <- grep("^p.", colnames(m.result.table))
rcols <- grep("^r.", colnames(m.result.table))
mcols <- grep("^m.", colnames(m.result.table))
m.result.table <- m.result.table[,c(pcols,rcols,mcols)]

pcols <- grep("^p.", colnames(p.result.table))
rcols <- grep("^r.", colnames(p.result.table))
mcols <- grep("^m.", colnames(p.result.table))
p.result.table <- p.result.table[,c(pcols,rcols,mcols)]

# Output ------------------------------------------------------------------
m.cols <- c("id", "symbol", "chr", "start", "end", "strand", "biotype")
m.output <- cbind(annot.mrna[,m.cols], m.result.table)
m.output$sig.mRNA_Age.Sex <- (m.output$p.mRNA_Age.Sex < 0.05)
m.output$sig.mRNA_Sex.Age <- (m.output$p.mRNA_Sex.Age < 0.05)

p.cols <- c("id", "gene_id", "symbol", "chr", "start", "end", "strand", "biotype")
p.output <- cbind(annot.protein[,p.cols], p.result.table)
p.output$sig.Prot_Age.Sex <- (p.output$p.Prot_Age.Sex < 0.05)
p.output$sig.Prot_Sex.Age <- (p.output$p.Prot_Sex.Age < 0.05)

dir <- "/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Anova_output/"
write.csv(m.output, paste0(dir, m.output.file), row.names=FALSE, quote = FALSE)
write.csv(p.output, paste0(dir, p.output.file), row.names=FALSE, quote = FALSE)

# How many are significant ?
dir <- "/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Anova_output/"
m.output <- read.csv(paste0(dir, m.output.file))
p.output <- read.csv(paste0(dir, p.output.file))

table(m.output$sig.mRNA_Age.Sex)
#> table(m.output$sig.mRNA_Age.Sex)
#
#FALSE  TRUE
#18859  3384
