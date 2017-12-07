library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggsci)
library(stringr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")

# Clean data
#df <- read.csv("./Phenotype/phenotypes/Lrp2RNAscope.csv")
#test <- t(as.data.frame(str_split(df$FileName_NucleiRaw, pattern = " - "))[2,])
#test <- t(as.data.frame(str_split(test, pattern = "-")))
#test[,2] <- str_pad(test[,2], 4, side = "left", pad = 0)
#df$sample <- paste0(test[,1],"-",test[,2])
#df$sample[df$sample == "DO-0954"] <- "DO-0952"  #fix typo
#df$sample[df$sample == "DO-1146"] <- "DO-1148"  #fix typo
#df$sample[df$sample == "D0-1263"] <- "DO-1263"  #fix typo
#df <- df[,c("sample", "Group", "Count_RNA", "Intensity_Protein", "Count_Nuceli", "FileName_NucleiRaw")]
#names(df) <- c("Sample", "Age", "RNAscope_Lrp2_count", "IF_LRP2_intensity", "Nuclei_count", "Filename")
#write.csv(df, "./Phenotype/phenotypes/Lrp2RNAscope_cleaned.csv", row.names = FALSE)

df <- read.csv("./Phenotype/phenotypes/Lrp2RNAscope_cleaned.csv")
df$ratio <- df$RNAscope_Lrp2_count / df$IF_LRP2_intensity

# Take log of genometic mean for each sample
animal <- as.data.frame(unique(df$Sample))
animal$Age <- NA
# Get age
for (i in 1:nrow(animal)){
  animal$Age[i] <- df[df$Sample == animal[,1][i],]$Age[1]
}
animal$log_geom <- NA
# Get geometic mean for each sample
for (i in 1:nrow(animal)){
  sub <- df[df$Sample == animal[,1][i],]
  animal$log_geom[i] <- log(prod(sub$ratio)^(1/nrow(sub)))
}
animal <- arrange(animal, Age)
names(animal)[1] <- "Sample"
animal$Age <- as.factor(as.character(animal$Age))
animal$Age <- factor(animal$Age,levels(animal$Age)[c(2,1)])

t.test(log_geom ~ Age, data = animal)$p.value
pdf("./Plot/Lrp2RNAscope_ttest.pdf", width = 8, height = 6)
ggplot(animal, aes(x = Age, y = log_geom, colour = Age)) +
      geom_boxplot() +
      theme_bw() +
      labs( title = "T-test: log(Geom.Mean(RNA count/protein)) at 6 and 18 months",
            subtitle = paste0("6 months n = 5, ", "\n",
                              "18 months n = 5, ", "\n",
                              "p-value = 0.5815"),
            y = "log(Geom.Mean(RNA count/protein))",
            x = "Age (months)") +
      guides( colour = FALSE) +
      scale_color_aaas()
dev.off()
