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

# Create function to do create the Jaccards distance df
makeJaccard <- function(file){

  # read file
  data <- read_csv(file)

  # create list of gene lists
  parsed_genelist <- vector("list", length = nrow(data))
  for(n in 1:nrow(data)){
    parsed_genelist[[n]] <- strsplit(data$geneID[n], "/")
  }

  # create empty df to fill
  Jaccard_df <- as.data.frame(matrix(data = NA, nrow = nrow(data), ncol = nrow(data)))
  rownames(Jaccard_df) <- colnames(Jaccard_df) <- data$ID

  # All vs All comparison
  for(x in 1:nrow(data)){
    print(x)
    for(y in 1:nrow(data)){
      AnB <- length(intersect(parsed_genelist[[x]][[1]], parsed_genelist[[y]][[1]]))
      AuB <- length(union(parsed_genelist[[x]][[1]], parsed_genelist[[y]][[1]]))
      Jaccard_df[x,y] <- AnB / AuB
    }
  }

  return(Jaccard_df)
}

# get list of all files, calcualte their Jaccard disance and make plots
list <- list.files("./results/Pathways", pattern ="QTLscan_chr", full.names = TRUE)
list <- list[grep(".csv$", list)]

# Loop through each file
for(i in 1:length(list)){

  # output base file name
  basename <- str_sub(list[i],,-5)
  print(i)
  print(basename)

  # calcualte Jaccard distance
  data <- read_csv(list[i])
  if(nrow(data) == 0){
    next
  }
  Jaccard_df <- makeJaccard(list[i])

  # plot heatmap
  if(nrow(Jaccard_df) > 1){
    set.seed(123)
    pdf(paste0(basename, "_heatmap.pdf"), width = 5, height = 5)
    heatmap.2(as.matrix(Jaccard_df),trace = "none", density.info = "none", key=FALSE, cexRow = 1, cexCol = 1, mar=c(6,6))
    dev.off()
  }

  # plot dendrogram
  if(nrow(Jaccard_df) > 2){
    hc <- hclust(dist(as.matrix(Jaccard_df)))
    pdf(paste0(basename, "_dendromap.pdf"), width = 5, height = 5)
    plot(hc)
    dev.off()
  }
}

# Cannot reduce, 0 or 1 pathway
