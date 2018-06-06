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

# Attempt 2 looking at correlations:
# annotate with snp marker names
table12 <- table12 %>% arrange(IntAgePos) %>% mutate( Marker = NA)
table12$Marker <- NA
for (i in 1:nrow(table12)){
  print(i)
  table12$Marker[i] <- snps[snps$chr == "12" & snps$bp %in% table12$IntAgePos[i],]$marker
}

# get all unique marker IDs
markers_unique <- unique(table12$Marker)

# subset genoprobs to only show markers_unique and get correlation
sub_geno <- genoprobs[,,markers_unique]

cor_data <- as.data.frame(matrix(data = NA, nrow= length(markers_unique), ncol = length(markers_unique)))
colnames(cor_data) <- markers_unique
rownames(cor_data) <- markers_unique

for(r in 1:length(markers_unique)){
  print(paste("Processing row: ",r))
  for(n in 1:length(markers_unique)){
    cor_data[r,n] <- cor(c(sub_geno[,,r]),c(sub_geno[,,n]), method ="pearson")
  }
}

library(reshape2)
tidy_cor <- melt(as.matrix(cor_data))
colnames(tidy_cor) <- c("Marker1", "Marker2", "PearsonCor")
TF <- c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)

#
png("./Plot/Chr12_markers_linkage.png", height = 1500, width = 1500)
ggplot(tidy_cor, aes(Marker1, Marker2, fill = PearsonCor))+
  geom_tile(colour = "white")+
  scale_fill_gradient2(low = "blue", high ="red", mid = "white", midpoint = 0.25, limit = c(min(tidy_cor$PearsonCor),max(tidy_cor$PearsonCor)), space = "Lab", name="Pearson\nCorrelation")+
  theme_minimal()+
  scale_x_discrete(breaks = markers_unique[rep(TF, times =53)])+
  scale_y_discrete(breaks = markers_unique[rep(TF, times =53)])+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
    size = 12, hjust = 1))+
  coord_fixed()
dev.off()

# Attempt 1 looking at alleles
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
