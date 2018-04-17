#qsub -v script=${R script name} Rsubmit_args.sh
# R/3.4.3
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney.RData")
options(dplyr.width = Inf)
#library(qtl2)
library(tidyverse)

# parameters
diffscan_dir <- "./SNPscan/diffscansnp_prot_addErk1_WBratio"
IntAge_output_file <- "./SNPscan/scan1snps_p_diffAgeInt_addErk1_WBratio_BestperGene.csv"

gatherdf <- function(dir, pattern){
  file_list <- list.files(dir,pattern,full.names=TRUE)

  #read and gather
  output <- NULL
  for(i in file_list){
    one <- read_csv(i) %>% mutate(chr = as.character(chr))
    output <- bind_rows(output,one)
  }

  # return df
  output <- output %>% arrange(id)
  return(output)
}

#gather
medLOD <- gatherdf(dir = diffscan_dir, pattern = "maxLODscan_batch_")
# write table
write_csv(medLOD, IntAge_output_file)

# get non mediated list
initialLOD <- readr::read_csv("./SNPscan/scan1snps_p_diffAgeInt_BestperGene_thr5.csv",
  guess_max = 4610) %>%
  filter(AgeIntChr == "7")

# check
if(!identical(medLOD$gene_id, initialLOD$gene_id)){
  initialLOD <- initialLOD[(initialLOD$id %in% medLOD$id),]
  medLOD <- medLOD[(medLOD$id %in% initialLOD$id),]
}

# combine and compare lod scores
compare <- initialLOD %>% mutate(
  Erk1_WBphos_MedChr = medLOD$AgeIntChr,
  Erk1_WBphos_MedPos = medLOD$AgeIntPos,
  Erk1_WBphos_MedLOD = medLOD$AgeIntLOD,
  LODdrop = AgeIntLOD - Erk1_WBphos_MedLOD
)

# write table
write_csv(compare, "./SNPscan/scan1snps_p_Erk1_WBratio_mediation_LODcomapre.csv")

# plot lod scores
pdf("./SNPscan/scan1snps_p_Erk1_WBratio_mediation_LODcompare.pdf", width = 9, height = 9)
ggplot(compare, aes(x=AgeIntLOD, y =Erk1_WBphos_MedLOD))+
  geom_point(alpha= 0.5)+
  geom_abline(intercept = 0, slope = 1, colour = "red")+
  geom_abline(intercept = -2, slope = 1, colour = "blue")+
  scale_x_continuous(name = "LOD score of Age Interactive proteome SNP scan", breaks = seq(0, 12, by = 1), labels = seq(0, 12, by = 1), limits = c(5,9))+
  scale_y_continuous(name = "LOD score of (X | p-Erk1:total-Erk1 ratio WB)", limits = c(2,9)) +
  theme_bw()+
  labs(title = "Protein SNPscan @Chr7 w/ p-Erk1:total-Erk1 ratio WB mediation",
       subtitle = paste0("Total genes: ", nrow(compare)))
dev.off()
