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
animal$log_geom_ratio <- NA
animal$log_geom_rna <- NA
animal$log_geom_protein <- NA
# Get geometic mean for each sample
for (i in 1:nrow(animal)){
  sub <- df[df$Sample == animal[,1][i],]
  animal$log_geom_ratio[i] <- log(prod(sub$ratio)^(1/nrow(sub)))
  animal$log_geom_rna[i] <- log(prod(sub$RNAscope_Lrp2_count)^(1/nrow(sub)))
  animal$log_geom_protein[i] <- log(prod(sub$IF_LRP2_intensity)^(1/nrow(sub)))
}
animal <- arrange(animal, Age)
names(animal)[1] <- "Sample"
animal$Age <- as.factor(as.character(animal$Age))
animal$Age <- factor(animal$Age,levels(animal$Age)[c(2,1)])

t.test(log_geom_ratio ~ Age, data = animal)$p.value
pdf("./Plot/Lrp2RNAscope_ttest.pdf", width = 8, height = 6)
ggplot(animal, aes(x = Age, y = log_geom_ratio, colour = Age)) +
      geom_boxplot() +
      theme_bw() +
      labs( title = "T-test: Megalin log(Geom.Mean(RNA count/protein)) at 6 and 18 months",
            subtitle = paste0("6 months n = 5, ", "\n",
                              "18 months n = 5, ", "\n",
                              "p-value = 0.5815"),
            y = "log(Geom.Mean(RNA count/protein))",
            x = "Age (months)") +
      guides( colour = FALSE) +
      scale_color_aaas()
dev.off()

pdf("./Plot/Lrp2RNAscope_mRNAvProt.pdf", width = 8, height = 6)
ggplot(animal, aes(x = log_geom_rna, y = log_geom_protein, colour = Age, fill = Age)) +
      geom_point(aes(shape = Age)) +
      labs( title = "Comparing Megalin RNA vs protein by Age",
            subtitle = paste0("6 months n = 5, ", "\n",
                              "18 months n = 5, "),
            y = "log(Geom.Mean(protein))",
            x = "log(Geom.Mean(RNA))") +
      scale_color_aaas()
dev.off()

# data redone------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(ggsci)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")

# get new data
df <- read.csv("./Phenotype/phenotypes/Lrp2RNAscope_trial2.csv")
df <- df[,c(1:6,11,12)]

# calculate
df$ProteinRNAratio <- df$Protein / df$CountRNA
df$ProteinRNAratio[df$ProteinRNAratio %in% c(Inf, NA)] <- NA
df$ProteinRNAratioPerTubArea <- df$ProteinRNAratio / df$Tubule.Area
df$ProteinRNAratioPerTubArea[df$ProteinRNAratioPerTubArea %in% c(Inf, NA)] <- NA
df$ProteinRNAratioPerNuclei <- df$ProteinRNAratio / df$X.Nuclei
df$ProteinRNAratioPerNuclei[df$ProteinRNAratioPerNuclei %in% c(Inf, NA)] <- NA
df$CountRNA[df$CountRNA %in% c(Inf, NA)] <- NA
df$Protein[df$Protein %in% c(Inf, NA)] <- NA

# Make new df and get age
animal <- as.data.frame(as.character(unique(df$Mouse)))
animal$Age <- NA
for (i in 1:nrow(animal)){
  animal$Age[i] <- as.character(df[df$Mouse == animal[,1][i],]$Group[1])
}
animal <- arrange(animal, Age)
names(animal)[1] <- "Sample"

# Calc log - geometric mean
animal$log_geom_ratio <- NA
animal$log_geom_ratioPerTub <- NA
animal$log_geom_ratioPerNuc <- NA
animal$log_geom_RNA <- NA
animal$log_geom_protein <- NA

for (i in 1:nrow(animal)){
  sub <- df[df$Mouse == animal[,1][i],]
  animal$log_geom_ratio[i] <- log(prod(sub$ProteinRNAratio, na.rm = TRUE)^(1/nrow(sub)))
  animal$log_geom_ratioPerTub[i] <- log(prod(sub$ProteinRNAratioPerTubArea, na.rm = TRUE)^(1/nrow(sub)))
  animal$log_geom_ratioPerNuc[i] <- log(prod(sub$ProteinRNAratioPerNuclei, na.rm = TRUE)^(1/nrow(sub)))
  animal$log_geom_RNA[i] <- log(prod(sub$CountRNA, na.rm = TRUE)^(1/nrow(sub)))
  animal$log_geom_protein[i] <- log(prod(sub$Protein, na.rm = TRUE)^(1/nrow(sub)))
}

# make factor and reorder levels
animal$Age <- as.factor(as.character(animal$Age))
animal$Age <- factor(animal$Age, levels(animal$Age)[c(2,1)])

# T-test
t.test(log_geom_ratio ~ Age, data = animal)$p.value
t.test(log_geom_ratioPerTub ~ Age, data = animal)$p.value
t.test(log_geom_ratioPerNuc ~ Age, data = animal)$p.value

animal$log_geom_RNA[animal$log_geom_RNA == -Inf] <- 0
pdf("./Plot/Lrp2RNAscope_mRNAvProt_try2_by_sample.pdf", width = 8, height = 6)
ggplot(animal, aes(x = log_geom_RNA, y = log_geom_protein, colour = Age, fill = Age)) +
      geom_point(aes(shape = Age)) +
      labs( title = "Comparing Megalin RNA vs protein by Age",
            subtitle = paste0("6 months n = 5, ", "\n",
                              "18 months n = 4, "),
            y = "log(Geom.Mean(protein))",
            x = "log(Geom.Mean(RNA))") +
      scale_color_aaas()
dev.off()

df$Group <- as.factor(as.character(df$Group))
pdf("./Plot/Lrp2RNAscope_mRNAvProt_try2_all.pdf", width = 8, height = 6)
ggplot(df, aes(x = log(CountRNA), y = log(Protein), colour = Group, fill = Group)) +
      geom_point(aes(shape = Group)) +
      labs( title = "Comparing Megalin RNA vs protein by Age",
            subtitle = paste0("6 months n = 206, ", "\n",
                              "18 months n = 81, "),
            y = "log(protein)",
            x = "log(RNA count)") +
      scale_color_aaas()
dev.off()
