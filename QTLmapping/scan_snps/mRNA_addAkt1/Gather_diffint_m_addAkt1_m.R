#qsub -v script=${R script name} Rsubmit_args.sh
# R/3.4.3
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney.RData")
options(dplyr.width = Inf)
#library(qtl2)
library(tidyverse)

# parameters
diffscan_dir <- "./SNPscan/diffscansnp_mrna_addAkt1_m"
IntAge_output_file <- "./SNPscan/scan1snps_m_diffAgeInt_addAkt1_m_BestperGene.csv"

gatherdf <- function(dir, pattern){
  dir <- diffscan_dir
  pattern <- "maxLOD"
  file_list <- list.files(dir,pattern,full.names=TRUE)

  #read and gather
  output <- NULL
  for(i in file_list){
    one <- readr::read_csv(i) %>%
      mutate(chr = as.character(chr),
        AgeIntChr = as.character(AgeIntChr))

    if(!is.null(output)){
      output$AgeIntChr <- as.character(output$AgeIntChr)
    }

    output <- bind_rows(output,one)
  }

  # return df
  output <- output %>% arrange(id)
  return(output)
}

diffLOD <- gatherdf(dir = diffscan_dir, pattern = "maxLODscan_batch_")

# write table
write_csv(diffLOD, IntAge_output_file)













initLOD <- readr::read_csv("./SNPscan/scan1snps_m_diffAgeInt_BestperGene_thr5.csv", guess_max = 4200) %>%
  filter(AgeIntChr == "12") %>% arrange(id)

compare <- initiLOD %>%





Chr12_list <- readr::read_csv("./SNPscan/scan1snps_m_diffAgeInt_BestperGene_thr5.csv", guess_max = 4200) %>%
  filter(AgeIntChr == "12")
sub_annot.mrna <- annot.mrna %>% filter(id %in% Chr12_list$id) #dim(sub_annot.mrna) 349 genes
sub_expr.mrna <- expr.mrna[,Chr12_list$id]



# output.file1 <- "./QTLscan/output/eQTLBestperGeneAddAkt1thr6_chr12.csv"
list_add <- read.csv(file = paste(output.file1), header = TRUE, stringsAsFactors = FALSE)
list_add <- arrange(list_add, id)

compare <- list[,colnames(list) %in% c("id", "symbol", "IntAgeChr", "IntAgePos", "IntAgeLODDiff")]
compare$addIntAgeChr <- list_add$IntAgeChr
compare$addIntAgePos <- list_add$IntAgePos
compare$addIntAgeLODDiff <- list_add$IntAgeLODDiff
compare <- compare[complete.cases(compare$addIntAgeChr),]
write.csv(compare, file="./QTLscan/output/eQTLintAkt1thr6_chr12.csv", row.names = FALSE)

# Plot Chr12 LOD scores
pdf("./QTLscan/output/plots/eQTL_Akt1Mediation_chr12_thr6.pdf", width = 9, heigh =9)
ggplot(compare, aes(x=IntAgeLODDiff,  y=addIntAgeLODDiff)) +
  geom_point(alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, color="red") +
  geom_abline(intercept = -2, slope = 1, color="blue") +
  scale_x_continuous(name = "LOD score Interactive age eQTL-diff", breaks = seq(0, 12, by = 1), labels = seq(0, 12, by = 1)) +
  scale_y_continuous(name = "LOD score (X | Akt1)", breaks = seq(0, 12, by = 1), labels = seq(0, 12, by = 1)) +
  theme_bw() +
  labs(title="eQTL Chr12 Genes Akt1 Mediation",
       subtitle = paste0("Chr 12 total: ", nrow(compare), " genes, threshold > 6 "))
dev.off()
