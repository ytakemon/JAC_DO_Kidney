# qsub -v script=QTLprot Rsubmit_args.sh

# R/3.4.1
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")
library(dplyr)

# parameters
addscan.dir <- "./QTLscan/addscan_prot/"
intscan.dir.Age <- "./QTLscan/intscan_prot_pbatch/Age/"
file.name <- function(i) paste0(annot.protein$id[i],"_",annot.protein$symbol[i],".rds")

LODtheshold_additive <- 7.5
LODthreshold_total <- 10.5
LODthreshold_diff <- 7.5

output.file1 <- "./QTLscan/output/pQTLBestperGene_nobatch.csv"
output.file2 <- "./QTLscan/output/pQTLAllAdditive_nobatch.csv"
output.file3 <- "./QTLscan/output/pQTLAllInteractiveAge_nobatch.csv"

output <- annot.protein[,c(1:6,10)]
output$AdditiveLOD <- output$AdditivePos <-  output$AdditiveChr <- NA
output$IntAgeLODFull <- output$IntAgeLODDiff <- output$IntAgePos <- output$IntAgeChr <- NA
output2 <- output3 <- list()

# for
for (i in 1:nrow(output)) {
  if (i %% 10 == 0) print(i)
  fit0 <- readRDS(paste0(addscan.dir,file.name(i)))
  fitAge <- readRDS(paste0(intscan.dir.Age, file.name(i)))

  # additive scan
  dt <- data.frame(AdditiveChr=sapply(strsplit(rownames(fit0),"_"), "[", 1),
                   AdditivePos=as.numeric(sapply(strsplit(rownames(fit0),"_"), "[", 2)),
                   AdditiveLOD=fit0[,1], stringsAsFactors=FALSE)
  dt2 <- dt %>% group_by(AdditiveChr) %>%
         summarize(AdditivePos = AdditivePos[which.max(AdditiveLOD)[1]],
                   AdditiveLOD = max(AdditiveLOD)) %>%
         arrange(-AdditiveLOD)

  output[i, c("AdditiveChr", "AdditivePos", "AdditiveLOD")] <- dt2[1,]
  to.be.added <-  as.data.frame(dt2 %>% filter(AdditiveLOD > LODtheshold_additive))
  if (nrow(to.be.added) > 0) output2[[i]] <- cbind(output[i,c(1:6)],to.be.added, row.names = NULL)

  # int. scan - Age
  stopifnot(rownames(fit0) == rownames(fitAge))
  dt <- data.frame(IntAgeChr=sapply(strsplit(rownames(fit0),"_"), "[", 1),
                   IntAgePos=as.numeric(sapply(strsplit(rownames(fit0),"_"), "[", 2)),
                   IntAgeLODDiff=fitAge[,1] - fit0[,1],
                   IntAgeLODFull = fitAge[,1],
                   stringsAsFactors=FALSE)
  dt2 <- dt %>% group_by(IntAgeChr) %>%
         summarize(IntAgePos = IntAgePos[which.max(IntAgeLODFull)[1]],
                   IntAgeLODDiff=IntAgeLODDiff[which.max(IntAgeLODFull)[1]],
                   IntAgeLODFull = max(IntAgeLODFull)) %>%
         arrange(-IntAgeLODDiff)
  output[i, c("IntAgeChr", "IntAgePos", "IntAgeLODDiff", "IntAgeLODFull")] <- dt2[1,]
  to.be.added <-  as.data.frame(dt2 %>% filter(IntAgeLODDiff > LODthreshold_diff, IntAgeLODFull > LODthreshold_total))
  if (nrow(to.be.added) > 0)  output3[[i]] <- cbind(output[i,c(1:6)],to.be.added, row.names = NULL)

}

# collect rows into one data frame
output2 <- do.call("rbind", output2)
output3 <- do.call("rbind", output3)

write.csv(output, file=output.file1, row.names=FALSE)
write.csv(output2, file=output.file2, row.names=FALSE)
write.csv(output3, file=output.file3, row.names=FALSE)
