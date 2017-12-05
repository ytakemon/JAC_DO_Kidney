library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggsci)
library(stringr)
library(reshape2)
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
hist(df$ratio) # skewed
hist(log(df$ratio)) # normalized
df$ratio_norm <- log(df$ratio)
# test
qqnorm(df$ratio_norm)
qqline(df$ratio_norm)
df$Age <- as.factor(df$Age)

fit <- lm(ratio_norm ~ Age + Sample, data = df)
summary(fit)

my_comparisons <- list( c(1, 2))
ratio_norm <- ggplot(df, aes(x  = Age, y = ratio_norm, colour = Age)) +
                geom_boxplot() +
                theme_bw() +
                labs( title = "log(RNA count/ protein) at 6 and 18 months",
                      subtitle = paste0("6 months n = 32, ", "\n",
                                        "18 months n = 31", "\n",
                                        "p-value = 1.564e-15"),
                      y = "log(RNA / protein)",
                      x = "Age (months)") +
                guides( colour = FALSE) +
                scale_color_aaas()

pdf("./Plot/Lrp2RNAscope_ttest.pdf", width = 8, height = 6)
ratio_norm
dev.off()

df$RNA_norm <- log(df$RNAscope_Lrp2_count)
RNA_norm <- ggplot(df, aes(x = Age, y = RNA_norm, colour = Age)) +
                  geom_boxplot() +
                  theme_bw() +
                  labs(title = "log(Lrp2 RNA) at 6 and 18 months",
                      subtitle = paste0("6 months n = 32, ", "\n",
                      "18 months n = 31"),
                      y = "log(RNA / protein)",
                      x = "Age (months)") +
                  stat_compare_means(label.y = 9.5, method = "t.test") +
                  guides( colour = FALSE) +
                  scale_color_aaas()

protein <- ggplot(df, aes(x = Age, y = IF_LRP2_intensity, colour = Age)) +
                  geom_boxplot() +
                  theme_bw() +
                  labs(title = "LRP2 protein at 6 and 18 months",
                      subtitle = paste0("6 months n = 32, ", "\n",
                      "18 months n = 31"),
                      y = "IF protein intensity",
                      x = "Age (months)") +
                  stat_compare_means(label.y = 0.03, method = "t.test") +
                  guides( colour = FALSE) +
                  scale_color_aaas()
