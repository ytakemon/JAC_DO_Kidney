# R/3.4.4
# Load each file and create a matrix of Jaccard distances:
# J(A,B) = |A n B| / |A u B|
#        = count will give % similarity

# load libraries
library(tidyverse)
library(reshape2)
library(gplots)

# set up directories
basedir <- "/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/"
setwd(basedir)

# Get list of files for comparison
slopeList <- list.files("./Pathways", pattern ="AgeSlope", full.names = TRUE)
slopeList <- slopeList[grep(glob2rx("*_clustered.csv$"), slopeList)]

anovaList <- list.files("./Pathways",pattern = "Age_clustered.csv", full.names = TRUE)

# check all 4 directions
for(direction in c("IncRNA_IncProt", "DecRNA_DecProt", "IncRNA_DecProt", "DecRNA_IncProt")){

  # Get sublist based on direction of change
  subList <- slopeList[grep(direction, slopeList)]

  # either mRNA or protein anova list
  for(type in c("mrna","protein")){

    # get subset of anova list for type
    if(type == "mrna"){
      alist <- anovaList[grep(glob2rx("*RNA*_clustered.csv$"), anovaList)]
    } else if(type == "protein"){
      alist <- anovaList[grep(glob2rx("*PROT*_clustered.csv$"), anovaList)]
    }

    # compare GO:bp
    if(length(subList[grep("GObp", subList)]) != 0){
      anova_GObp <- read.csv(alist[grep("GObp", alist)], stringsAsFactors = FALSE)
      slope_GObp <- read.csv(subList[grep("GObp", subList)], stringsAsFactors = FALSE)
      GObp_intersect <- intersect(anova_GObp$ID, slope_GObp$ID)
      slope_GObp$Overlap_ANOVA <- FALSE

      if(length(GObp_intersect) != 0){
        for(term in GObp_intersect){
          slope_GObp[slope_GObp$ID == term, ]$Overlap_ANOVA <- TRUE
        }
      }
      basename_bp <- str_sub(subList[grep("GObp", subList)],,-5)
      write.csv(slope_GObp, paste0(basename_bp,"_compared.csv"), row.names = FALSE)
    }

    # compare GO:cc
    if(length(subList[grep("GOcc", subList)]) != 0){
      anova_GOcc <- read.csv(alist[grep("GOcc", alist)], stringsAsFactors = FALSE)
      slope_GOcc <- read.csv(subList[grep("GOcc", subList)], stringsAsFactors = FALSE)
      GOcc_intersect <- intersect(anova_GOcc$ID, slope_GOcc$ID)
      slope_GOcc$Overlap_ANOVA <- FALSE

      if(length(GOcc_intersect) != 0){
        for(term in GOcc_intersect){
          slope_GOcc[slope_GOcc$ID == term, ]$Overlap_ANOVA <- TRUE
        }
      }
      basename_cc <- str_sub(subList[grep("GOcc", subList)],,-5)
      write.csv(slope_GOcc, paste0(basename_cc,"_compared.csv"), row.names = FALSE)
    }

    # compare KEGG
    if(length(subList[grep("KEGG", subList)]) != 0){
      anova_KEGG <- read.csv(alist[grep("KEGG", alist)], stringsAsFactors = FALSE)
      slope_KEGG <- read.csv(subList[grep("KEGG", subList)], stringsAsFactors = FALSE)
      KEGG_intersect <- intersect(anova_KEGG$ID, slope_KEGG$ID)
      slope_KEGG$Overlap_ANOVA <- FALSE

      if(length(KEGG_intersect) != 0){
        for(term in KEGG_intersect){
          slope_KEGG[slope_KEGG$ID == term, ]$Overlap_ANOVA <- TRUE
        }
      }
      basename_KEGG <- str_sub(subList[grep("KEGG", subList)],,-5)
      write.csv(slope_KEGG, paste0(basename_KEGG,"_compared.csv"), row.names = FALSE)
    }
  }
}
