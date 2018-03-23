#qsub -v script=${R script name} Rsubmit_args.sh
# R/3.4.3
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney.RData")
library(dplyr)

# parameters
addscan.dir <- "./QTLscan/addscan_prot_perAge/6mo/"
file.name <- function(i) paste0(annot.protein$id[i],"_",annot.protein$symbol[i],".rds")

output.file <- "./QTLscan/output/6mo_pQTLBestperGene.csv"

output <- annot.protein[,c(1:6,10)]
output$AdditiveLOD <- output$AdditivePos <-  output$AdditiveChr <- NA

# for
for (i in 1:nrow(output)) {
  if (i %% 10 == 0) print(i)
  fit0 <- readRDS(paste0(addscan.dir,file.name(i)))

  # additive scan
  dt <- data.frame(AdditiveChr=sapply(strsplit(rownames(fit0),"_"), "[", 1),
                   AdditivePos=as.numeric(sapply(strsplit(rownames(fit0),"_"), "[", 2)),
                   AdditiveLOD=fit0[,1], stringsAsFactors=FALSE)
  dt2 <- dt %>% group_by(AdditiveChr) %>%
         summarize(AdditivePos = AdditivePos[which.max(AdditiveLOD)[1]],
                   AdditiveLOD = max(AdditiveLOD)) %>%
         arrange(-AdditiveLOD)

  output[i, c("AdditiveChr", "AdditivePos", "AdditiveLOD")] <- dt2[1,]
}
write.csv(output, file=output.file, row.names=FALSE)
