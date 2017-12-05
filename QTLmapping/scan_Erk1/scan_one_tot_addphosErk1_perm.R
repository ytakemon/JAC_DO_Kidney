# qsub -v script=scan_one_tot_addphosErk1_perm Rsubmit_args.sh

library(qtl2)
library(qtl2convert)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO1045_kidney.Rdata")
erk1 <- read.delim("./Phenotype/phenotypes/JAC_WB_kidney_ERK.txt")
pheno <- "Total_ERK1"

# Cleanup data and subset to match samples
# match pheno to samples
erk1 <- erk1[complete.cases(erk1[,pheno]) & complete.cases(erk1[,"Phospho_ERK1"]),]
erk1 <- erk1[erk1$ID %in% samples$Mouse.ID, ]
erk1$duplicated <- (duplicated(erk1$ID) | duplicated(erk1$ID, fromLast = TRUE))
erk1 <- erk1[erk1$duplicated == FALSE,]
rownames(erk1) <- erk1$ID

# match samples to pheno
sub_samples <- samples[rownames(samples) %in% rownames(erk1),]

#match genoprobs to samples
sub_genoprobs <- genoprobs[rownames(genoprobs) %in% rownames(sub_samples),,]

# Check identical
identical(rownames(erk1), rownames(sub_samples))
identical(rownames(erk1), rownames(sub_genoprobs))

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(sub_genoprobs, MM_snps, pos_column = "pos")
K <- calc_kinship(probs, "loco", cores=1)
MM_snps$chr <- as.character(MM_snps$chr)
MM_snps$chr[MM_snps$chr=="20"] <- "X"
snps <- MM_snps[dimnames(sub_genoprobs)[[3]],]
map <- map_df_to_list(map = snps, pos_column = "pos")

sub_samples$add <- erk1$Phospho_ERK1
addcovar <- model.matrix(~ Sex + Generation + Cohort.Age.mo + add, data = sub_samples)

# read lod
file_name <- paste0("./QTLscan/addscan_phenotype/", pheno, "addPhos.rds")
lod <- readRDS(file=file_name)

# permutation test
perm <- scan1perm(genoprobs = probs,
                  pheno = erk1[,pheno, drop = FALSE],
                  kinship = K,
                  addcovar = addcovar[,-1],
                  n_perm = 1000,
                  cores = 10,
                  reml = TRUE)

# save perms
file_name <- paste0("./QTLscan/addscan_phenotype/", pheno, "addPhos_perm.rds")
perm <- saveRDS(file_name)

# plot
#plot(lod, map)
#title(main = "TITLE")
#abline(h = 4, col = "red")
#abline(h = 5)
