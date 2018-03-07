library(tidyverse)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")
output.file1 <- "./QTLscan/output/eQTLBestperGeneAddAkt1thr6_chr12.csv"
list <- read.csv("./QTLscan/output/Threshold6_eQTL_intAge.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == 12, ]
list <- arrange(list, id)

# output.file1 <- "./QTLscan/output/eQTLBestperGeneAddAkt1thr6_chr12.csv"
list_add <- read.csv(file = paste(output.file1), header = TRUE, stringsAsFactors = FALSE)
list_add <- arrange(list_add, id)

compare <- list[,colnames(list) %in% c("id", "symbol", "IntAgeChr", "IntAgePos", "IntAgeLODDiff")]
compare$addIntAgeChr <- list_add$IntAgeChr
compare$addIntAgePos <- list_add$IntAgePos
compare$addIntAgeLODDiff <- list_add$IntAgeLODDiff
compare <- compare[complete.cases(compare$addIntAgeChr),]

ggplot(compare, aes(x = addIntAgePos, y = addIntAgeLODDiff)) +
  geom_point() +
  geom_abline(intercept = 6, slope = 0, color="blue")

ggplot(compare, aes(x = IntAgePos, y = IntAgeLODDiff)) +
  geom_point() +
  geom_abline(intercept = 6, slope = 0, color="blue")


output.file1 <- paste0("./QTLscan/output/pQTLBestperGene_ERK1Ratio_thr6_chr7.csv")
list <- read.csv("./QTLscan/output/Threshold6_pQTL_intAge_pbatch.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == 7, ]
list <- arrange(list, id)

list_add <- read.csv(file = paste(output.file1), header = TRUE, stringsAsFactors = FALSE)
list_add <- arrange(list_add, id)

compare2 <- list[,colnames(list) %in% c("id", "gene_id", "symbol", "IntAgeChr", "IntAgePos", "IntAgeLODDiff")]
compare2$addIntAgeChr <- list_add$IntAgeChr
compare2$addIntAgePos <- list_add$IntAgePos
compare2$addIntAgeLODDiff <- list_add$IntAgeLODDiff
compare2 <- compare2[complete.cases(compare2$addIntAgeChr),]

ggplot(compare2, aes(x = addIntAgePos, y = addIntAgeLODDiff)) +
  geom_point() +
  geom_abline(intercept = 6, slope = 0, color="blue")

ggplot(compare2, aes(x = IntAgePos, y = IntAgeLODDiff)) +
  geom_point() +
  geom_abline(intercept = 6, slope = 0, color="blue")
