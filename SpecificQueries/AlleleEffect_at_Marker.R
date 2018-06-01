library(qtl2)
library(qtl2convert)
library(tidyverse)
#library(ggsci)
options(dplyr.width = Inf) #override column limit
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")

# get table to use
IntAge_output_file <- "./QTLscan/scanBestMarker_mrna/BestMarker_BestperGene_mrna.csv"
table <- read_csv(IntAge_output_file)
table12 <- table %>% filter(IntAgeChr == 12, IntAgeLODDiff > 7)
ggplot(table12, aes(x= IntAgePos, y= IntAgeLODDiff)) +
geom_point()

# Tabulate frequency of markers that appear on chr12
freq <- table(table12$IntAgePos)
max(freq) # max is 10
length(freq) # There are 325
freq <- freq[freq >5]
length(freq) # There are 12
markerpos <- as.numeric(names(freq))

# Get allele distribution at these clusters
GetAlleleDist <- function(chr, pos){

  # get marker id
  id <- snps[snps$chr == chr & snps$bp %in% pos,]$marker

  # create df
  df <- data.frame(Mouse.ID = annot.samples$Mouse.ID)

  for (m in 1:length(id)){
    alleles <- as.data.frame(genoprobs[,,id[m]])

    Geno <- NULL
    for (i in 1:nrow(alleles)){
      # Find Hets
      if(length(alleles[i,][alleles[i,] > 0.3]) == 2){
        hets <- names(alleles[i,][which(alleles[i,] > 0.3)])
        geno <- paste0(hets[1],hets[2])
      } else if(length(alleles[i,][alleles[i,] > 0.3]) == 1){
        hom <- names(alleles[i,][which(alleles[i,] > 0.3)])
        geno <- paste0(hom,hom)
      } else{
        geno <- "ZZ"
      }
      Geno <- c(Geno,geno)
    }

    df[,paste(id[m])] <- Geno
  }
  return(df)
}

id <- snps[snps$chr == 12 & snps$bp %in% markerpos,]$marker
GetAlleles <- GetAlleleDist(12,markerpos)
Collect <- GetAlleles %>% gather(MarkerPOS, Allele, colnames(GetAlleles)[2:ncol(GetAlleles)])
Collect$Allele <- as.factor(Collect$Allele)
Collect$MarkerPOS <- as.factor(Collect$MarkerPOS)
Collect$MarkerPOS <- factor(Collect$MarkerPOS, levels = id)



pdf("./Plot/Chr12_clustermarkers_allele.pdf", height = 4, width = 8)
ggplot(Collect, aes(Allele, fill = MarkerPOS, group = MarkerPOS))+
  geom_bar()+
  labs(title= "Chr12 marker freq > 5")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
