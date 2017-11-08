# R/3.4.1
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")
library(ggplot2)
library(dplyr)
library(scales)

# Get list of genes with trans eQTL
list <- read.csv("./QTLscan/output/Threshold6_eQTL_intAge.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == 12, ]

# parameters
addscan.dir <- "./QTLscan/addscan_mrna_WB_PanTotAkt1/"
intscan.dir.Age <- "./QTLscan/intscan_mrna_WB_PanTotAkt1/"
output.file1 <- "./QTLscan/output/eQTLBestperGeneAddWB_PanTotAkt1thr6.csv"

annot.mrna <- annot.mrna[annot.mrna$id %in% list$id,]
output <- annot.mrna[,c(1:5,9)]
output$AdditiveLOD <- output$AdditivePos <-  output$AdditiveChr <- NA
output$IntAgeLODFull <- output$IntAgeLODDiff <- output$IntAgePos <- output$IntAgeChr <- NA

file.name <- function(i) paste0(output$id[i],"_",output$symbol[i],".rds")

# for
for (i in 1:nrow(output)) {
  if (i %% 10 == 0) print(i)

  if(file.exists(paste0(addscan.dir,file.name(i))) && file.exists(paste0(intscan.dir.Age, file.name(i)))){
    fit0 <- readRDS(paste0(addscan.dir,file.name(i)))
    fitAge <- readRDS(paste0(intscan.dir.Age, file.name(i)))
  }else{
    next
  }

  # additive scan
  dt <- data.frame(AdditiveChr=sapply(strsplit(rownames(fit0),"_"), "[", 1),
                   AdditivePos=as.numeric(sapply(strsplit(rownames(fit0),"_"), "[", 2)),
                   AdditiveLOD=fit0[,1], stringsAsFactors=FALSE)
  dt2 <- dt %>% group_by(AdditiveChr) %>%
         summarize(AdditivePos = AdditivePos[which.max(AdditiveLOD)[1]],
                   AdditiveLOD = max(AdditiveLOD)) %>%
         arrange(-AdditiveLOD)

  output[i, c("AdditiveChr", "AdditivePos", "AdditiveLOD")] <- dt2[1,]

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
}

# collect rows into one data frame
write.csv(output, file=output.file1, row.names=FALSE)

# compare:
list <- read.csv("./QTLscan/output/Threshold6_eQTL_intAge.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == 12, ]
list <- arrange(list, id)

list_add <- read.csv("./QTLscan/output/eQTLBestperGeneAddWB_PanTotAkt1thr6", header = TRUE, stringsAsFactors = FALSE)
list_add <- arrange(list_add, id)

compare <- list[,colnames(list) %in% c("id", "symbol", "IntAgeChr", "IntAgeLODDiff")]
compare$addIntAgeChr <- list_add$IntAgeChr
compare$addIntAgeLODDiff <- list_add$IntAgeLODDiff
compare$change <- !(compare$IntAgeChr == compare$addIntAgeChr)
compare <- compare[complete.cases(compare$addIntAgeChr),]
write.csv(compare, file="./QTLscan/output/eQTLintPanTotAkt1thr6.csv")

# Plot LOD score
pdf("./QTLscan/output/plots/eQTL_PanTotAkt1Mediation_thr6.pdf", width = 9, heigh =9)
ggplot(compare, aes(x=IntAgeLODDiff,  y=addIntAgeLODDiff, colour = change)) +
  geom_point(alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, color="red") +
  guides(colour=guide_legend(title = "Mediation")) +
  xlab("LOD score Interactive age eQTL-diff") +
  ylab("LOD score (X | Pan-total-AKT1)") +
  theme_bw() +
  labs(title="Pan-total-AKT1 eQTL Chr12 Genes Mediation",
       subtitle = paste0("Chr 12 total: ", nrow(compare), " genes, mediated: ", table(compare$change)[[2]], " genes, threshold > 6 "))
dev.off()
