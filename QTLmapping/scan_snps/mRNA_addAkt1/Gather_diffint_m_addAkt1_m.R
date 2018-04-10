#qsub -v script=${R script name} Rsubmit_args.sh
# R/3.4.3
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney.RData")
options(dplyr.width = Inf)
#library(qtl2)
library(tidyverse)

# parameters
diffscan_dir <- "./SNPscan/diffscansnp_mrna_addAkt1_m"
gatherdf <- function(dir, pattern){
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

# gather all scans
medLOD <- gatherdf(dir = diffscan_dir, pattern = "maxLODscan_batch_")

# get non mediated list
initialLOD <- readr::read_csv("./SNPscan/scan1snps_m_diffAgeInt_BestperGene_thr5.csv", guess_max = 4200) %>%
  filter(AgeIntChr == "12") %>% arrange(id)

# check
identical(medLOD$id, initialLOD$id)

# combine and compare lod scores
compare <- initialLOD %>% mutate(
  Akt1MedChr = medLOD$AgeIntChr,
  Akt1MedPos = medLOD$AgeIntPos,
  Akt1MedLOD = medLOD$AgeIntLOD,
  LODdrop = AgeIntLOD - Akt1MedLOD
)

# write table
write_csv(compare, "./SNPscan/scan1snps_m_Akt1_m_mediation_LODcomapre.csv")

# plot lod scores
pdf("./SNPscan/scan1snps_m_Akt1_m_mediation_LODcompare.pdf", width = 9, height = 9)
ggplot(compare, aes(x=AgeIntLOD, Akt1MedLOD))+
  geom_point(alpha= 0.5)+
  geom_abline(intercept = 0, slope = 1, colour = "red")+
  geom_abline(intercept = -2, slope = 1, colour = "blue")+
  scale_x_continuous(name = "LOD score of Age Interactive transcriptome SNP scan", breaks = seq(0, 12, by = 1), labels = seq(0, 12, by = 1))+
  scale_y_continuous(name = "LOD score of (X | Akt1 mRNA)") +
  theme_bw()+
  labs(title = "mRNA SNPscan @Chr12 w/ Akt1 mRNA mediation",
       subtitle = paste0("Total genes: ", nrow(compare)))
dev.off()
