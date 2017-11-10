# R/3.4.1
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")
library(ggplot2)
library(dplyr)
library(scales)

# Get list of genes with trans pQTL
list <- read.csv("./QTLscan/output/Threshold6_pQTL_intAge.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == 15, ]

# parameters
addscan.dir <- "./QTLscan/addscan_prot_batch/"
intscan.dir.Age <- "./QTLscan/intscan_prot_batch/"
output.file1 <- "./QTLscan/output/pQTLBestperGeneAddBatchthr6_chr15.csv"

annot.protein <- annot.protein[annot.protein$id %in% list$id,]
output <- annot.protein[,c(1:6,10)]
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

  # find max lod score of whole genome
  # dt2 <- dt %>% group_by(AdditiveChr) %>%
  #       summarize(AdditivePos = AdditivePos[which.max(AdditiveLOD)[1]],
  #                  AdditiveLOD = max(AdditiveLOD)) %>%
  #        arrange(-AdditiveLOD)
  # output[i, c("AdditiveChr", "AdditivePos", "AdditiveLOD")] <- dt2[1,]

  # find max lod score of chr 15 ONLY!
  dt2 <- dt[dt$AdditiveChr == "15",]
  dt2 <- dt %>% group_by(AdditiveChr) %>%
        summarize(AdditivePos = AdditivePos[which.max(AdditiveLOD)[1]],
                   AdditiveLOD = max(AdditiveLOD)) %>%
         arrange(-AdditiveLOD)
  output[i, c("AdditiveChr", "AdditivePos", "AdditiveLOD")] <- dt2[dt2$AdditiveChr == "15",]

  # int. scan - Age
  stopifnot(rownames(fit0) == rownames(fitAge))
  dt <- data.frame(IntAgeChr=sapply(strsplit(rownames(fit0),"_"), "[", 1),
                   IntAgePos=as.numeric(sapply(strsplit(rownames(fit0),"_"), "[", 2)),
                   IntAgeLODDiff=fitAge[,1] - fit0[,1],
                   IntAgeLODFull = fitAge[,1],
                   stringsAsFactors=FALSE)

  # find max lod score of whole genome
  # dt2 <- dt %>% group_by(IntAgeChr) %>%
  #        summarize(IntAgePos = IntAgePos[which.max(IntAgeLODFull)[1]],
  #                  IntAgeLODDiff=IntAgeLODDiff[which.max(IntAgeLODFull)[1]],
  #                  IntAgeLODFull = max(IntAgeLODFull)) %>%
  #        arrange(-IntAgeLODDiff)
  # output[i, c("IntAgeChr", "IntAgePos", "IntAgeLODDiff", "IntAgeLODFull")] <- dt2[1,]

  # find max lod score of chr15 ONLY!
  dt2 <- dt %>% group_by(IntAgeChr) %>%
         summarize(IntAgePos = IntAgePos[which.max(IntAgeLODFull)[1]],
                   IntAgeLODDiff=IntAgeLODDiff[which.max(IntAgeLODFull)[1]],
                   IntAgeLODFull = max(IntAgeLODFull)) %>%
         arrange(-IntAgeLODDiff)
  output[i, c("IntAgeChr", "IntAgePos", "IntAgeLODDiff", "IntAgeLODFull")] <- dt2[dt2$IntAgeChr == "15",]
}

# collect rows into one data frame
write.csv(output, file=output.file1, row.names=FALSE)

# compare:
list <- read.csv("./QTLscan/output/Threshold6_pQTL_intAge.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == 15, ]
list <- arrange(list, id)

# output.file1 <- "./QTLscan/output/pQTLBestperGeneAddAkt1thr6_chr15.csv"
list_add <- read.csv(file = paste(output.file1), header = TRUE, stringsAsFactors = FALSE)
list_add <- arrange(list_add, id)

if (identical(list$id, list_add$id)){
  compare <- list[,colnames(list) %in% c("id", "gene_id", "symbol", "IntAgeChr", "IntAgePos", "IntAgeLODDiff")]
  compare$addIntAgeChr <- list_add$IntAgeChr
  compare$addIntAgePos <- list_add$IntAgePos
  compare$addIntAgeLODDiff <- list_add$IntAgeLODDiff
  compare <- compare[complete.cases(compare$addIntAgeChr),]
  write.csv(compare, file="./QTLscan/output/pQTLintBatchthr6_chr15.csv", row.names = FALSE)
} else {
  print("Lists do not match")
}


# Plot Chr15 LOD scores
pdf("./QTLscan/output/plots/pQTL_BatchbMediation_chr15_thr6.pdf", width = 9, heigh =9)
ggplot(compare, aes(x=IntAgeLODDiff,  y=addIntAgeLODDiff)) +
  geom_point(alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, color="red") +
  geom_abline(intercept = -2, slope = 1, color="blue") +
  scale_x_continuous(name = "LOD score Interactive age pQTL-diff", breaks = seq(0, 15, by = 1), labels = seq(0, 15, by = 1)) +
  scale_y_continuous(name = "LOD score (X | Proteomics Batch & Channel)", breaks = seq(0, 12, by = 1), labels = seq(0, 12, by = 1)) +
  theme_bw() +
  labs(title="pQTL Chr15 Genes Proteomics Batch & Channel Mediation",
       subtitle = paste0("Chr 15 total: ", nrow(compare), " genes, threshold > 6 "))
dev.off()

# Plot LOD score *retired*
# pdf("./QTLscan/output/plots/pQTL_Akt1Mediation_thr6.pdf", width = 9, heigh =9)
# ggplot(compare, aes(x=IntAgeLODDiff,  y=addIntAgeLODDiff, colour = change)) +
#   geom_point(alpha=0.5) +
#   geom_abline(intercept = 0, slope = 1, color="red") +
#   guides(colour=guide_legend(title = "Mediation")) +
#   xlab("LOD score Interactive age pQTL-diff") +
#   ylab("LOD score (X | Akt1)") +
#   theme_bw() +
#   labs(title="Akt1 pQTL Chr12 Genes Mediation",
#        subtitle = paste0("Chr 12 total: ", nrow(compare), " genes, mediated: ", table(compare$change)[[2]], " genes, threshold > 6 "))
# dev.off()
