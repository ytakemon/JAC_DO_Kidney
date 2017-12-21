library(ggplot2)
library(dplyr)
library(ggsci)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
dropDock1 <- read.csv("./QTLscan/output/pERK1_upstream/protein/pERK1_upstream_prot_mediator_ENSMUSP00000081531_Dock1_compare_chr7_thr6.csv")
dropGlrx3 <- read.csv("./QTLscan/output/pERK1_upstream/protein/pERK1_upstream_prot_mediator_ENSMUSP00000066621_Glrx3_compare_chr7_thr6.csv")
dropPhosERK1 <- read.csv("./QTLscan/output/pQTLint_ERK1Phos_chr7_thr6.csv")

compare <- dropPhosERK1[,c(3,6,9)]
compare$gene <- "PhosERK1"
temp <- dropDock1[,c(3,6,9)]
temp$gene <- "Dock1"
compare <- rbind(compare, temp)

pdf(paste0("./QTLscan/output/plots/pERK1_mediation_upstream_compare_Dock1_change.pdf"), width = 20, heigh =18)
ggplot(compare, aes(x = IntAgeLODDiff, y = addIntAgeLODDiff, group = symbol, colour = gene)) +
  geom_point(alpha=0.5) +
  geom_line(alpha=0.2, colour = "black")+
  geom_abline(intercept = 0, slope = 1, color="red") +
  geom_abline(intercept = -2, slope = 1, color="blue") +
  scale_x_continuous( name = "LOD score Interactive age pQTL-diff",
                      breaks = seq(0, 15, by = 1),
                      labels = seq(0, 15, by = 1)) +
  scale_y_continuous( name = paste0("LOD score (X | mediator)"),
                      breaks = seq(0, 15, by = 1),
                      labels = seq(0, 15, by = 1)) +
  theme_bw() +
  scale_color_aaas() +
  labs(title=paste0("pERK1 mediation comp. upstream Dock1 mediation"),
       subtitle = paste0(""))
dev.off()

compare <- dropPhosERK1[,c(3,6,9)]
compare$gene <- "PhosERK1"
temp <- dropGlrx3[,c(3,6,9)]
temp$gene <- "Glrx3"
compare <- rbind(compare, temp)

pdf(paste0("./QTLscan/output/plots/pERK1_mediation_upstream_compare_Glrx3_change.pdf"), width = 20, heigh =18)
ggplot(compare, aes(x = IntAgeLODDiff, y = addIntAgeLODDiff, group = symbol, colour = gene)) +
  geom_point(alpha=0.5) +
  geom_line(alpha=0.2, colour = "black")+
  geom_abline(intercept = 0, slope = 1, color="red") +
  geom_abline(intercept = -2, slope = 1, color="blue") +
  scale_x_continuous( name = "LOD score Interactive age pQTL-diff",
                      breaks = seq(0, 15, by = 1),
                      labels = seq(0, 15, by = 1)) +
  scale_y_continuous( name = paste0("LOD score (X | mediator)"),
                      breaks = seq(0, 15, by = 1),
                      labels = seq(0, 15, by = 1)) +
  theme_bw() +
  scale_color_aaas() +
  labs(title=paste0("pERK1 mediation comp. upstream Glrx3 mediation"),
       subtitle = paste0(""))
dev.off()

identical(dropPhosERK1$id, dropDock1$id, dropGlrx3$id)
compare <- dropPhosERK1
compare$Dock1 <- dropDock1$addIntAgeLODDiff
compare$Glrx3 <- dropGlrx3$addIntAgeLODDiff
write.csv(compare, "./QTLscan/output/pERK1_upstream/pERK1_Dock1_Glrx3_directcomp.csv", row.names = FALSE)

pdf(paste0("./QTLscan/output/plots/pERK1_mediation_upstream_compare_Dock1_direct.pdf"), width = 20, heigh =18)
ggplot(compare, aes(x = addIntAgeLODDiff, y = Dock1)) +
  geom_point(alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, color="red") +
  geom_abline(intercept = -2, slope = 1, color="blue") +
  geom_abline(intercept = 2, slope = 1, color="blue") +
  scale_x_continuous( name = "LOD score (X | pERK1)",
                      breaks = seq(0, 15, by = 1),
                      labels = seq(0, 15, by = 1)) +
  scale_y_continuous( name = paste0("LOD score (X | Dock1)"),
                      breaks = seq(0, 15, by = 1),
                      labels = seq(0, 15, by = 1)) +
  theme_bw() +
  labs(title=paste0("pERK1 mediation comp. upstream Dock1 mediation"),
       subtitle = paste0(""))
dev.off()

pdf(paste0("./QTLscan/output/plots/pERK1_mediation_upstream_compare_Glrx3_direct.pdf"), width = 20, heigh =18)
ggplot(compare, aes(x = addIntAgeLODDiff, y = Glrx3)) +
  geom_point(alpha=0.5) +
  geom_abline(intercept = 0, slope = 1, color="red") +
  geom_abline(intercept = -2, slope = 1, color="blue") +
  geom_abline(intercept = 2, slope = 1, color="blue") +
  scale_x_continuous( name = "LOD score (X | pERK1)",
                      breaks = seq(0, 15, by = 1),
                      labels = seq(0, 15, by = 1)) +
  scale_y_continuous( name = paste0("LOD score (X | Glrx3)"),
                      breaks = seq(0, 15, by = 1),
                      labels = seq(0, 15, by = 1)) +
  theme_bw() +
  labs(title=paste0("pERK1 mediation comp. upstream Glrx3 mediation"),
       subtitle = paste0(""))
dev.off()

compare$Dock1Diff <- compare$addIntAgeLODDiff - compare$Dock1
compare$Glrx3Diff <- compare$addIntAgeLODDiff - compare$Glrx3
write.csv(compare, "./QTLscan/output/pERK1_upstream/pERK1_Dock1_Glrx3_change2.csv", row.names = FALSE)

pdf(paste0("./QTLscan/output/plots/pERK1_mediation_upstream_compare_Dock1_change2.pdf"), width = 20, heigh =18)
ggplot(compare, aes(x = addIntAgeLODDiff, y = Dock1Diff)) +
  geom_point(alpha=0.5) +
  geom_abline(intercept = 0, slope = 0, color="red") +
  geom_abline(intercept = -2, slope = 0, color="blue") +
  geom_abline(intercept = 2, slope = 0, color="blue") +
  scale_x_continuous( name = "LOD score (X | pERK1)",
                      breaks = seq(-15, 15, by = 1),
                      labels = seq(-15, 15, by = 1)) +
  scale_y_continuous( name = paste0("Change by Dock1 mediation"),
                      breaks = seq(-15, 15, by = 1),
                      labels = seq(-15, 15, by = 1)) +
  theme_bw() +
  labs(title=paste0("pERK1 mediation comp. upstream Dock1 mediation"),
       subtitle = paste0("pERK1 LOD - Dock1 LOD"))
dev.off()

pdf(paste0("./QTLscan/output/plots/pERK1_mediation_upstream_compare_Glrx3_change2.pdf"), width = 20, heigh =18)
ggplot(compare, aes(x = addIntAgeLODDiff, y = Glrx3Diff)) +
  geom_point(alpha=0.5) +
  geom_abline(intercept = 0, slope = 0, color="red") +
  geom_abline(intercept = -2, slope = 0, color="blue") +
  geom_abline(intercept = 2, slope = 0, color="blue") +
  scale_x_continuous( name = "LOD score (X | pERK1)",
                      breaks = seq(-15, 15, by = 1),
                      labels = seq(-15, 15, by = 1)) +
  scale_y_continuous( name = paste0("Change by Glrx3 mediation"),
                      breaks = seq(-15, 15, by = 1),
                      labels = seq(-15, 15, by = 1)) +
  theme_bw() +
  labs(title=paste0("pERK1 mediation comp. upstream Glrx3 mediation"),
       subtitle = paste0("pERK1 LOD - Glrx3 LOD"))
dev.off()

# Venn Diagram-----------------------------------------------------
dropDock1 <- read.csv("./QTLscan/output/pERK1_upstream/protein/pERK1_upstream_prot_mediator_ENSMUSP00000081531_Dock1_compare_chr7_thr6.csv")
dropGlrx3 <- read.csv("./QTLscan/output/pERK1_upstream/protein/pERK1_upstream_prot_mediator_ENSMUSP00000066621_Glrx3_compare_chr7_thr6.csv")
dropPhosERK1 <- read.csv("./QTLscan/output/pQTLint_ERK1Phos_chr7_thr6.csv")

dropDock1$diff <- dropDock1$IntAgeLODDiff - dropDock1$addIntAgeLODDiff
dropDock1$lod2 <- dropDock1$diff >= 2
dropGlrx3$diff <- dropGlrx3$IntAgeLODDiff - dropGlrx3$addIntAgeLODDiff
dropGlrx3$lod2 <- dropGlrx3$diff >= 2
dropPhosERK1$diff <- dropPhosERK1$IntAgeLODDiff - dropPhosERK1$addIntAgeLODDiff
dropPhosERK1$lod2 <- dropPhosERK1$diff >= 2

if (identical(dropPhosERK1$id, dropDock1$id, dropGlrx3$id)){
  compare <- dropPhosERK1
  compare$drop2Dock1 <- dropDock1$lod2
  compare$drop2Glrx3 <- dropGlrx3$lod2
  compare$drop2pERK_Dock1 <- compare$lod2 & compare$drop2Dock1 & TRUE
  compare$drop2pERK_Glrx3 <- compare$lod2 & compare$drop2Glrx3 & TRUE
  compare$drop2Dock1_Glrx3 <- compare$drop2Dock1 & compare$drop2Glrx3
  compare$drop2trio <- compare$lod2 & compare$drop2Dock1 & compare$drop2Glrx3 & TRUE
} else (
  print("can't compare, check id")
)



table(compare$lod2)
table(compare$drop2Dock1)
table(compare$drop2Glrx3)
table(compare$drop2pERK_Dock1)
table(compare$drop2pERK_Glrx3)
table(compare$drop2Dock1_Glrx3)
table(compare$drop2trio)
