#qsub -v script=${R script name} Rsubmit_args.sh
# R/3.4.3
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney.RData")
options(dplyr.width = Inf)
#library(qtl2)
library(tidyverse)

# parameters
diffscan_dir <- "./QTLscan/scanBestMarker_protein/maxLODscan_Erk1_WBratio"
gatherdf <- function(dir, pattern){
  file_list <- list.files(dir,pattern,full.names=TRUE)

  #read and gather
  output <- NULL
  for(i in file_list){
    one <- readr::read_csv(i) %>%
      mutate(chr = as.character(chr),
        IntAgeChr = as.character(IntAgeChr))

    if(!is.null(output)){
      output$IntAgeChr <- as.character(output$IntAgeChr)
    }

    output <- bind_rows(output,one)
  }

  # return df
  output <- output %>% arrange(id)
  return(output)
}

# gather all scans
medLOD <- gatherdf(dir = diffscan_dir, pattern = "maxLODscan_batch_")

# get non mediated list
initialLOD <- readr::read_csv("./QTLscan/scanBestMarker_protein/BestMarker_BestperGene_protein_thr8.csv", guess_max = 4200) %>% filter(IntAgeChr == "7") %>% arrange(id)

# check
identical(medLOD$id, initialLOD$id)

# combine and compare lod scores
compare <- initialLOD %>% mutate(
  Erk1MedChr = medLOD$IntAgeChr,
  Erk1MedPos = medLOD$IntAgePos,
  Erk1MedLOD = medLOD$IntAgeLODDiff,
  LODdrop = IntAgeLODDiff - Erk1MedLOD
)

# write table
write_csv(compare, "./QTLscan/scanBestMarker_protein/maxLODscan_Erk1_WBratio/scanBestMarker_p_Erk1_WBratio_mediation_LODcomapre.csv")

# plot lod scores
pdf("./QTLscan/scanBestMarker_protein/maxLODscan_Erk1_WBratio/scan1snps_p_Erk1_WBratio_mediation_LODcompare.pdf", width = 9, height = 9)
ggplot(compare, aes(x=IntAgeLODDiff, Erk1MedLOD))+
  geom_point(alpha= 0.5)+
  geom_abline(intercept = 0, slope = 1, colour = "red")+
  geom_abline(intercept = -2, slope = 1, colour = "blue")+
  scale_x_continuous(name = "LOD score of Age Interactive proteome scan", breaks = seq(0, 12, by = 1), labels = seq(0, 12, by = 1))+
  scale_y_continuous(name = "LOD score of (X | pERK1/totalErk1)") +
  theme_bw()+
  labs(title = "Protein QTLscan @Chr7 w/ pErk1/totalErk1 protein mediation",
       subtitle = paste0("Total genes: ", nrow(compare)))
dev.off()
