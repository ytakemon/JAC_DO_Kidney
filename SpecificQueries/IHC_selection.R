library(ggplot2)
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")

# Selection quality:
Sex <- "M"
Age <- c(6, 18)
Num <- 3 # Number of animals to select

# 
