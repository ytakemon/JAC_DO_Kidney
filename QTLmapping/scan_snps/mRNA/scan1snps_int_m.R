plist <- as.numeric(commandArgs(trailingOnly = TRUE))

cat(paste("Mapping", length(plist), "genes", "\n"))
print(Sys.time())
print(plist)

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
plist <- plist[plist<=ncol(expr.mrna)]

# create output file for lod score harvest
output <- annot.mrna[plist,-ncol(annot.mrna)]
output$FullLOD  <- output$FullPos <- output$FullChr <- NA

for (p in plist) {
  # print message
  cat("Scanning ",which(p==plist)," out of ",length(plist),"\n")
  file_name <- paste0("./SNPscan/intscansnp_mrna/", annot.mrna$id[p], "_", annot.mrna$symbol[p], ".rds")

  addcovar <- model.matrix(~ Sex + Age + Generation, data=annot.samples)
  intcovar <- model.matrix(~ Age, data=annot.samples)
  # Perform scan1snps
  # If query_func is given, but start and end are empty, it should calcualte for all chromosomes.
  snpsOut <- scan1snps(genoprobs=probs,
               map = map,
               kinship=Glist,
               pheno=expr.mrna[,p],
               addcovar=addcovar[,-1],
               intcovar=intcovar[,-1],
               query_func = query_func,
               chr =c(1:19,"X"),
               start = 0,
               end = 200,
               keep_all_snps = FALSE,
               cores=20, reml=TRUE)

  # save highest lod object
  rsid <-rownames(snpsOut$lod)[which.max(snpsOut$lod)]
  # assign to output file
  output[which(p==plist),]$FullLOD <- snpsOut$lod[which.max(snpsOut$lod)]
  output[which(p==plist),]$FullChr <- snpsOut$snpinfo[snpsOut$snpinfo$snp_id == rsid,]$chr
  output[which(p==plist),]$FullPos <- snpsOut$snpinfo[snpsOut$snpinfo$snp_id == rsid,]$pos
}

write.csv(output, file = paste0("./SNPscan/intscansnp_mrna/maxLODscan_batch_",plist[1],".csv"),row.names = FALSE)
print(Sys.time())

# Following warning messages will appear and its fine:
#Warning messages:
#1: In scan1snps(genoprobs = probs, map = map, kinship = Glist, pheno = expr.mrna[,  :
#  If length(chr) > 1, start end end are ignored.
