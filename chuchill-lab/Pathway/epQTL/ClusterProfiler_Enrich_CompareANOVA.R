# R/3.4.4
# Load each file and create a matrix of Jaccard distances:
# J(A,B) = |A n B| / |A u B|
#        = count will give % similarity

# load libraries
library(tidyverse)
library(reshape2)
library(gplots)

# set up directories
basedir <- "/projects/churchill-lab/projects/JAC/Takemon_DO_crosssectional_kidney/"
setwd(basedir)

# Get list of files for comparison
QTLlist <- list.files("./results/Pathways", pattern ="QTLscan", full.names = TRUE)

anovaList <- list.files("./results/Pathways",pattern = "Age_clustered.csv", full.names = TRUE)

# check all both QTLtypes
for(QTLtype in c("mrna","protein")){

  # Get sublist based on QTLtype of change
  subList <- QTLlist[grep(QTLtype, QTLlist)]

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
      qlt_GObp <- read.csv(subList[grep("GObp", subList)], stringsAsFactors = FALSE)
      # if exists
      if(nrow(qlt_GObp) != 0){
        GObp_intersect <- intersect(anova_GObp$ID, qlt_GObp$ID)
        qlt_GObp$Overlap_ANOVA <- FALSE
        # if any overlap
        if(length(GObp_intersect) != 0){
          for(term in GObp_intersect){
            qlt_GObp[qlt_GObp$ID == term, ]$Overlap_ANOVA <- TRUE
          }
        }
        # save
        basename_bp <- str_sub(subList[grep("GObp", subList)],,-5)
        if(type == "mrna"){
          write.csv(qlt_GObp, paste0(basename_bp,"_mrnaANOVA_compared.csv"), row.names = FALSE)
        } else if(type == "protein"){
          write.csv(qlt_GObp, paste0(basename_bp,"_proteinANOVA_compared.csv"), row.names = FALSE)
        }

      }
    }

    # compare GO:cc
    if(length(subList[grep("GOcc", subList)]) != 0){
      anova_GOcc <- read.csv(alist[grep("GOcc", alist)], stringsAsFactors = FALSE)
      qlt_GOcc <- read.csv(subList[grep("GOcc", subList)], stringsAsFactors = FALSE)
      # if exists
      if(nrow(qlt_GOcc) != 0){
        GOcc_intersect <- intersect(anova_GOcc$ID, qlt_GOcc$ID)
        qlt_GOcc$Overlap_ANOVA <- FALSE
        # if any overlap
        if(length(GOcc_intersect) != 0){
          for(term in GOcc_intersect){
            qlt_GOcc[qlt_GOcc$ID == term, ]$Overlap_ANOVA <- TRUE
          }
        }
        # save
        basename_cc <- str_sub(subList[grep("GOcc", subList)],,-5)
        if(type == "mrna"){
          write.csv(qlt_GOcc, paste0(basename_cc,"_mrnaANOVA_compared.csv"), row.names = FALSE)
        } else if(type == "protein"){
          write.csv(qlt_GOcc, paste0(basename_cc,"_proteinANOVA_compared.csv"), row.names = FALSE)
        }
      }
    }

    # compare KEGG
    if(length(subList[grep("KEGG", subList)]) != 0){
      anova_KEGG <- read.csv(alist[grep("KEGG", alist)], stringsAsFactors = FALSE)
      qlt_KEGG <- read.csv(subList[grep("KEGG", subList)], stringsAsFactors = FALSE)
      # if exists
      if(nrow(qlt_KEGG) != 0){
        KEGG_intersect <- intersect(anova_KEGG$ID, qlt_KEGG$ID)
        qlt_KEGG$Overlap_ANOVA <- FALSE
        # if any overlap
        if(length(KEGG_intersect) != 0){
          for(term in KEGG_intersect){
            qlt_KEGG[qlt_KEGG$ID == term, ]$Overlap_ANOVA <- TRUE
          }
        }
        # save
        basename_KEGG <- str_sub(subList[grep("KEGG", subList)],,-5)
        if(type == "mrna"){
          write.csv(qlt_KEGG, paste0(basename_KEGG,"_mrnaANOVA_compared.csv"), row.names = FALSE)
        } else if (type == "protein"){
          write.csv(qlt_KEGG, paste0(basename_KEGG,"_proteinANOVA_compared.csv"), row.names = FALSE)
        }
      }
    }
  }
}
