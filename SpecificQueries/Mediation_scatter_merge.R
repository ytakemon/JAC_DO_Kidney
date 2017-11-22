library(dplyr)
library(ggplot2)

setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("RNAseq_data/DO188b_kidney_noprobs.RData")

# Get slope data
file <- "kidney_anova_slope_output.csv"
data <- read.csv(paste0("./Anova_output/",file[[1]]), header = T)
# For only slopes and sigcols
data <- data[, c(1:9,11,41,43)]

# Get mediation data
MAP <- read.csv("./QTLscan/output/pQTLint_ERK1Ratio_chr7_thr6.csv", stringsAsFactors = FALSE)
Akt1 <- read.csv("./QTLscan/output/eQTLintAkt1thr6_chr12.csv", stringsAsFactors = FALSE)
# Get lod drop > 2
MAP$drop <- MAP$IntAgeLODDiff - MAP$addIntAgeLODDiff
Akt1$drop <- Akt1$IntAgeLODDiff - Akt1$addIntAgeLODDiff
MAP$drop2 <- (MAP$drop >= 2)
Akt1$drop2 <- (Akt1$drop >= 2)
# subset drop 2 genes only
MAP <- MAP[MAP$drop2 == TRUE,]
Akt1 <- Akt1[Akt1$drop2 == TRUE,]

# Integrate into total data
data$MAP_drop2 <- data$gene_id %in% MAP$gene_id
data$Akt1_drop2 <- data$gene_id %in% Akt1$id

# sigonly
data_sig <- data[data$p.mRNA_Age.Sex <= 0.05,]
data_sig <- data_sig[data_sig$p.Prot_Age.Sex <= 0.05,]
write.csv(data_sig, file = "./SpecificQ/Slop_medGene.csv")

# MAP total
pdf("./Plot/Total_slope_ERK1ratio_medGenes.pdf", width = 8, heigh =6)
ggplot() +
      geom_point( data = data,
                  aes(x = m.mRNA_Age.Sex,
                      y = m.Prot_Age.Sex,
                      colour = MAP_drop2),
                  alpha = 0.7) +
      scale_colour_manual( values = c("#C2C0BE", "#F86A08")) +
      guides( colour = guide_legend(title = "pERK1/ERK1 ratio Mediated Genes")) +
      scale_x_continuous("mRNA-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      scale_y_continuous("Protein-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs( title = "pERK1/ERK1 ratio Mediated genes in total list view",
            subtitle = "LOD drop > 2") +
      theme(panel.background = element_blank(),
            panel.border = element_rect( colour = "black", fill = NA))
dev.off()


#MAP sig
pdf("./Plot/Sig_slope_ERK1ratio_medGenes.pdf", width = 8, heigh =6)
ggplot() +
      geom_point( data = data_sig,
                  aes(x = m.mRNA_Age.Sex,
                      y = m.Prot_Age.Sex,
                      colour = MAP_drop2),
                  alpha = 0.7) +
      geom_text(data = data_sig,
                aes(x = m.mRNA_Age.Sex,
                    y = m.Prot_Age.Sex,
                    label = ifelse(data_sig$MAP_drop2 == TRUE, as.character(data_sig$symbol),""),
                    hjust=0,
                    vjust=0))+
      scale_colour_manual( values = c("#C2C0BE", "#F86A08")) +
      guides( colour = guide_legend(title = "pERK1/ERK1 ratio Mediated Genes")) +
      scale_x_continuous("mRNA-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      scale_y_continuous("Protein-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs( title = "pERK1/ERK1 ratio Mediated genes in significant list view",
            subtitle = "LOD drop > 2") +
      theme(panel.background = element_blank(),
            panel.border = element_rect( colour = "black", fill = NA))
dev.off()

# Akt total
pdf("./Plot/Total_slope_Akt1_medGenes.pdf", width = 8, heigh =6)
ggplot() +
      geom_point( data = data,
                  aes(x = m.mRNA_Age.Sex,
                      y = m.Prot_Age.Sex,
                      colour = Akt1_drop2),
                  alpha = 0.7) +
      scale_colour_manual( values = c("#C2C0BE", "#F86A08")) +
      guides( colour = guide_legend(title = "Akt1 Mediated Genes")) +
      scale_x_continuous("mRNA-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      scale_y_continuous("Protein-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs( title = "Akt1 Mediated genes in total list view",
            subtitle = "LOD drop > 2") +
      theme(panel.background = element_blank(),
            panel.border = element_rect( colour = "black", fill = NA))
dev.off()

# Akt sig
pdf("./Plot/Sig_slope_Akt1_medGenes.pdf", width = 8, heigh =6)
ggplot() +
      geom_point( data = data_sig,
                  aes(x = m.mRNA_Age.Sex,
                      y = m.Prot_Age.Sex,
                      colour = Akt1_drop2),
                  alpha = 0.7) +
      geom_text(data = data_sig,
                aes(x = m.mRNA_Age.Sex,
                    y = m.Prot_Age.Sex,
                    label = ifelse(data_sig$Akt1_drop2 == TRUE, as.character(data_sig$symbol),""),
                    hjust=0,
                    vjust=0))+
      scale_colour_manual( values = c("#C2C0BE", "#F86A08")) +
      guides( colour = guide_legend(title = "Akt1 Mediated Genes")) +
      scale_x_continuous("mRNA-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      scale_y_continuous("Protein-Age.Sex",
                         breaks = seq(-1, 1, by = 0.1)) +
      geom_vline(xintercept = 0) +
      geom_hline(yintercept = 0) +
      labs( title = "Akt1 Mediated genes in significant list view",
            subtitle = "LOD drop > 2") +
      theme(panel.background = element_blank(),
            panel.border = element_rect( colour = "black", fill = NA))
dev.off()


# Test drop --------------------------------------------------------------------
# MAPK3
chr <- "7"
g <- "ERK1Ratio"
# parameters
addscan.dir <- paste0("./QTLscan/addscan_prot_", g ,"/")
intscan.dir.Age <-  paste0("./QTLscan/intscan_prot_", g, "/")
output.file1 <- paste0("./QTLscan/output/pQTLBestperGene_", g, "_thr6_chr7.csv")
list <- read.csv("./QTLscan/output/Threshold6_pQTL_intAge_pbatch.csv", header = TRUE, stringsAsFactors = FALSE)
list <- list[list$IntAgeChr == chr, ]
list <- arrange(list, id)
list_add <- read.csv(file = paste(output.file1), header = TRUE, stringsAsFactors = FALSE)
list_add <- arrange(list_add, id)

if (identical(list$id, list_add$id)){
  compare <- list[,colnames(list) %in% c("id", "gene_id", "symbol", "IntAgeChr", "IntAgePos", "IntAgeLODDiff")]
  compare$addIntAgeChr <- list_add$IntAgeChr
  compare$addIntAgePos <- list_add$IntAgePos
  compare$addIntAgeLODDiff <- list_add$IntAgeLODDiff
  compare <- compare[complete.cases(compare$addIntAgeChr),]
  write.csv(compare, file=paste0("./QTLscan/output/pQTLint_", g, "_chr7_thr6.csv"), row.names = FALSE)
} else {
  print("Lists do not match")
}

compare$MAP_drop2 <- compare$gene_id %in% MAP$gene_id
# Plot Chr15 LOD scores
ggplot(compare, aes(x=IntAgeLODDiff,  y=addIntAgeLODDiff, colour = MAP_drop2)) +
  geom_point(alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, color="red") +
  geom_abline(intercept = -2, slope = 1, color="blue") +
  scale_x_continuous( name = "LOD score Interactive age pQTL-diff",
                      breaks = seq(0, 15, by = 1),
                      labels = seq(0, 15, by = 1)) +
  scale_y_continuous( name = paste0("LOD score (X | d mRNA)"),
                      breaks = seq(0, 15, by = 1),
                      labels = seq(0, 15, by = 1)) +
  theme_bw()
