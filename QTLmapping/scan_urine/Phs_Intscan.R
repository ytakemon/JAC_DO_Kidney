# qsub -v script=Phs_Intscan Rsubmit_args.sh
library(qtl2convert)
library(qtl2)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO1045_kidney.Rdata")

# Subset pheontype: phs
pheno <- Upheno[Upheno$study == "Cross-sectional",]
pheno <- pheno[-1,]
pheno_cr <- pheno[,c("Mouse.ID", "cr.u.6", "cr.u.12", "cr.u.18")]
pheno_phs <- pheno[,c("Mouse.ID", "phs.u.6", "phs.u.12", "phs.u.18")]

# colapse to one column and combine
pheno_cr <- pheno_cr %>% mutate(
              cr.u.all = coalesce(cr.u.6, cr.u.12, cr.u.18)
            )
pheno_phs <- pheno_phs %>% mutate(
              phs.u.all = coalesce(phs.u.6, phs.u.12, phs.u.18)
            )

pheno <- pheno_cr[,c("Mouse.ID","cr.u.all")]
pheno$phs.u.all <- pheno_phs$phs.u.all
pheno <- pheno[!is.na(pheno[,2] & pheno[,3]),]
rownames(pheno) <- pheno$Mouse.ID

# Subset dataset
genoprobs <- genoprobs[pheno$Mouse.ID,,]
samples <- samples[pheno$Mouse.ID,]
samples$cr.u.all <- pheno$cr.u.all

# Identified missing age in samples
samples[which(is.na(samples$Cohort.Age.mo)),]
#"DO-0835" NA
#"DO-0906" NA
#"DO-0930" 18
#"DO-0943" 18
#"DO-0947" 18
#"DO-1044" NA
#"DO-1189" 12
#"DO-1246" 6
#"DO-1317" 6
samples[which(is.na(samples$Cohort.Age.mo)),]$Cohort.Age.mo <- c(NA, NA, 18, 18, 18, NA, 12, 6, 6)
samples <- samples[-which(is.na(samples$Cohort.Age.mo)),]

# resubset: Total should be 446
genoprobs <- genoprobs[as.character(samples$Mouse.ID),,]
pheno <- pheno[as.character(samples$Mouse.ID),]

# prepare data for qtl2
MM_snps$chr <- as.character(MM_snps$chr)
probs <- probs_doqtl_to_qtl2(genoprobs, MM_snps, pos_column = "pos")
K <- calc_kinship(probs, type = "loco", cores = 20)
map <- map_df_to_list(map = MM_snps, pos_column = "pos")
samples$cr <- pheno[,2]
addcovar <- model.matrix(~ Sex + Cohort.Age.mo + Generation + Cohort + cr.u.all , data = samples)
intcovar <- model.matrix(~ Cohort.Age.mo, data = samples)

# scan
lod <- scan1(genoprobs=probs,
             kinship=K,
             pheno=as.data.frame(pheno$phs.u.all, row.names = rownames(pheno)),
             addcovar=addcovar[,-1],
             intcovar=intcovar[,-1],
             cores=20,
             reml=TRUE)

# save lod
saveRDS(lod, file = "./QTLscan/addscan_urine/Intscan_phs_all.rds")

perm <- scan1perm(genoprobs=probs,
                     kinship=K,
                     pheno=as.data.frame(pheno$phs.u.all, row.names = rownames(pheno)),
                     addcovar=addcovar[,-1],
                     intcovar=intcovar[,-1],
                     cores=20,
                     n_perm = 1000,
                     reml = TRUE)

# save permutation
saveRDS(perm, file = "./QTLscan/addscan_urine/Intperm_phs_all.rds")

# Get coef
# get max lod
chr <- max(lod, map)$chr
# calc coef
coef <- scan1coef(genoprobs = probs[,chr],
                  kinship = K[chr],
                  pheno = as.data.frame(pheno$phs.u.all, row.names = rownames(pheno)),
                  addcovar = addcovar[,-1],
                  intcovar=intcovar[,-1],
                  reml = TRUE)
# save coef
saveRDS(coef, file = "./QTLscan/addscan_urine/Intcoef_phs_all.rds")

# Get genes in lod peak interval
query_variants <- create_variant_query_func("./qtl2_sqlite/cc_variants.sqlite")
peak_Mbp <- max(lod, map)$pos
peak_chr <- max(lod, map)$chr

out_snps <- scan1snps(genoprobs = probs,
                      map = map,
                      pheno = as.data.frame(pheno$phs.u.all, row.names = rownames(pheno)),
                      kinship =K[[peak_chr]],
                      addcovar = addcovar[,-1],
                      intcovar=intcovar[,-1],
                      query_func=query_variants,
                      chr=peak_chr,
                      start=peak_Mbp-1,
                      end=peak_Mbp+1,
                      keep_all_snps=TRUE,
                      cores = 20)
saveRDS(out_snps, file = "./QTLscan/addscan_urine/Intsnps_phs_all.rds")
