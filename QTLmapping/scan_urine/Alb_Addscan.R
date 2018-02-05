# qsub -v script=Alb_Addscan Rsubmit_args.sh
library(qtl2convert)
library(qtl2)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO1045_kidney.Rdata")

# Subset pheontype: Alb
pheno <- Upheno[Upheno$study == "Cross-sectional",]
pheno <- pheno[-1,]
pheno_cr <- pheno[,c("Mouse.ID", "cr.u.6", "cr.u.12", "cr.u.18")]
pheno_ma <- pheno[,c("Mouse.ID", "ma.u.6", "ma.u.12", "ma.u.18")]

# colapse to one column and combine
pheno_cr <- pheno_cr %>% mutate(
              cr.u.all = coalesce(cr.u.6, cr.u.12, cr.u.18)
            )
pheno_ma <- pheno_ma %>% mutate(
              ma.u.all = coalesce(ma.u.6, ma.u.12, ma.u.18)
            )

pheno <- pheno_cr[,c("Mouse.ID","cr.u.all")]
pheno$ma.u.all <- pheno_ma$ma.u.all
pheno <- pheno[!is.na(pheno[,2] & pheno[,3]),]
rownames(pheno) <- pheno$Mouse.ID

# Subset dataset
genoprobs <- genoprobs[pheno$Mouse.ID,,]
samples <- samples[pheno$Mouse.ID,]
samples$cr.u.all <- pheno$cr.u.all

# Identified missing age in samples
samples[which(is.na(samples$Cohort.Age.mo)),]
#DO-0906 NA
#DO-0930 18mo.
#DO-0943 18mo.
#DO-1189 12mo.
samples[which(is.na(samples$Cohort.Age.mo)),]$Cohort.Age.mo <- c(NA, 18, 18, 12)
samples <- samples[-which(is.na(samples$Cohort.Age.mo)),]

# resubset: Total should be 308
genoprobs <- genoprobs[as.character(samples$Mouse.ID),,]
pheno <- pheno[as.character(samples$Mouse.ID),]

# prepare data for qtl2
MM_snps$chr <- as.character(MM_snps$chr)
probs <- probs_doqtl_to_qtl2(genoprobs, MM_snps, pos_column = "pos")
K <- calc_kinship(probs, type = "loco", cores = 20)
map <- map_df_to_list(map = MM_snps, pos_column = "pos")
samples$cr <- pheno[,2]
addcovar <- model.matrix(~ Sex + Cohort.Age.mo + Generation + Cohort + cr.u.all , data = samples)

# scan
lod <- scan1(genoprobs=probs,
             kinship=K,
             pheno=as.data.frame(pheno$ma.u.all, row.names = rownames(pheno)),
             addcovar=addcovar[,-1],
             cores=20,
             reml=TRUE)
# save lod
saveRDS(lod, file = "./QTLscan/addscan_urine/Addscan_alb_all.rds")

# Get perm
perm <- scan1perm(genoprobs=probs,
                     kinship=K,
                     pheno=as.data.frame(pheno$ma.u.all, row.names = rownames(pheno)),
                     addcovar=addcovar[,-1],
                     cores=20,
                     n_perm = 1000,
                     reml = TRUE)
# save permutation
saveRDS(perm, file = "./QTLscan/addscan_urine/Addperm_alb_all.rds")

# Get coef
# get max lod
chr <- max(lod, map)$chr
pos <- max(lod, map)$pos
# calc coef
coef <- scan1coef(genoprobs = probs[,chr],
                  kinship = K[chr],
                  pheno = as.data.frame(pheno$ma.u.all, row.names = rownames(pheno)),
                  addcovar = addcovar[,-1],
                  reml = TRUE)
# save coef
saveRDS(coef, file = "./QTLscan/addscan_urine/Addcoef_alb_all.rds")

# Get genes in lod peak interval

query_variants <- create_variant_query_func("./qtl2_sqlite/cc_variants.sqlite")













###
scansnp <- scan1snps( genoprobs = probs,
                      kinship = K,
                      pheno = as.data.frame(pheno$ma.u.all, row.names = rownames(pheno)),
                      map = map,
                      addcovar = addcovar[,-1],
                      query_func = create_variant_query_func(2, pos - 1, pos + 1),
                      chr = chr,
                      start= pos - 1,
                      end = pos + 1,
                      keep_all_snps = TRUE,
                      reml = TRUE)

scansnp
