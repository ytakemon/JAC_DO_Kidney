plist <- as.numeric(commandArgs(trailingOnly = TRUE))

cat(paste("Mapping", length(plist), "genes", "\n"))
print(Sys.time())
print(plist)

library(tidyverse)
library(qtl2)
library(qtl2convert)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")
snpdb <- "/hpcdata/gac/resource/CCsnps/cc_variants_v2.sqlite"

# prepare data for qtl2
probs <- probs_doqtl_to_qtl2(genoprobs, snps, pos_column = "bp")
map <- map_df_to_list(map = snps, pos_column = "bp")
query_func <- create_variant_query_func(snpdb)

# check list to not exceed list
plist <- plist[plist<=ncol(expr.protein)]

# create output file for lod score harvest
output <- annot.protein[plist,]
output$AgeIntLOD  <- output$AgeIntPos <- output$AgeIntChr <- NA

for (p in plist) {
  # print message
  cat("Scanning ",which(p==plist)," out of ",length(plist),"\n")

  addcovar <- model.matrix(~ Sex + Age + Generation + Protein.Batch + Protein.Channel, data=annot.samples)
  intcovar <- model.matrix(~ Age, data=annot.samples)
  # Perform scan1snps
  # If query_func is given, but start and end are empty, it should calcualte for all chromosomes.
  snpsOut_add <- scan1snps(genoprobs=probs,
                    map = map,
                    kinship=Glist,
                    pheno=expr.protein[,p],
                    addcovar=addcovar[,-1],
                    query_func = query_func,
                    chr =c(1:19,"X"),
                    start = 0,
                    end = 200,
                    keep_all_snps = FALSE,
                    cores=20, reml=TRUE)

  snpsOut_full <- scan1snps(genoprobs=probs,
                    map = map,
                    kinship=Glist,
                    pheno=expr.protein[,p],
                    addcovar=addcovar[,-1],
                    intcovar=intcovar[,-1],
                    query_func = query_func,
                    chr =c(1:19,"X"),
                    start = 0,
                    end = 200,
                    keep_all_snps = FALSE,
                    cores=20, reml=TRUE)

  # line up both snp scan
  if(!identical(rownames(snpsOut_add$lod), rownames(snpsOut_full$lod)) & nrow(snpsOut_add$lod) == nrow(snpsOut_full$lod)){
    stop("either numer of rows or rsids do not match!")
  }

  # find max lod diff
  LODcomp <- as.data.frame(snpsOut_add$lod) %>% rename(
    add_lod = pheno1
  ) %>% mutate(
    rsid = rownames(snpsOut_add$lod),
    full_lod = snpsOut_full$lod[,1],
    diff_lod = full_lod - add_lod,
  ) %>% filter(diff_lod == max(diff_lod, na.rm = TRUE))

  # assign to output file
  output[which(p==plist),]$AgeIntLOD <- LODcomp$diff_lod[1]
  output[which(p==plist),]$AgeIntChr <- snpsOut_add$snpinfo[snpsOut_add$snpinfo$snp_id == LODcomp$rsid[1],]$chr
  output[which(p==plist),]$AgeIntPos <- snpsOut_add$snpinfo[snpsOut_add$snpinfo$snp_id == LODcomp$rsid[1],]$pos
}

write.csv(output, file = paste0("./SNPscan/diffscansnp_prot/maxLODscan_batch_",plist[1],".csv"),row.names = FALSE)
print(Sys.time())

# Following warning messages will appear and its fine:
#Warning messages:
#1: In scan1snps(genoprobs = probs, map = map, kinship = Glist, pheno = expr.protein[,  :
#  If length(chr) > 1, start end end are ignored.


# problem found: p
[1] 5112


# plist  <- was okay
 [1] 5170 5171 5172 5173 5174 5175 5176 5177 5178 5179 5180 5181 5182 5183 5184
[16] 5185 5186 5187 5188 5189 5190 5191 5192 5193 5194 5195 5196 5197 5198 5199
[31] 5200
