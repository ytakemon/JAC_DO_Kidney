# R/3.4.1
library(ggplot2)
library(dplyr)
library(grid)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
eQTL_best <- read.csv("./QTLscan/output/eQTLBestperGene.csv")

# Using eQTL_best to create plot
# need to reorder "chr" factors
chr_full <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X","Y", "MT")
eQTL_best$chr <- factor(eQTL_best$chr, levels = chr_full)
eQTL_best$IntAgeChr <- factor(eQTL_best$IntAgeChr, levels= chr_full)
eQTL_best$IntSexChr <- factor(eQTL_best$IntSexChr, levels= chr_full)
# Subset out chr 1-19,X from data
chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X")
eQTL_best <- eQTL_best[eQTL_best$chr %in% chr, ]

# Combine chr&positions into one numerical verctor
eQTL_best$Gene_point <- paste(eQTL_best$chr, eQTL_best$start, sep = ".")
eQTL_best$IntAge_point <- paste(eQTL_best$IntAgeChr, eQTL_best$IntAgePos, sep = ".")
eQTL_best$IngSex_point <- paste(eQTL_best$IntSexChr, eQTL_best$IntSexPos, sep = ".")

# Set LOD threshold
LODthreshold_diff <- 8

# Plot Interactive-Age eQTLs
# Subset Int-Age LOD above total/full and diff
Int_age <- eQTL_best[(eQTL_best$IntAgeLODDiff > LODthreshold_diff),] # above diff threshold

# Annotate Interactive-Age postion with genes and save file for sharing
save_int_age <- Int_age[,c("id", "symbol","chr","start","end", "biotype", "IntAgeChr","IntAgePos","IntAgeLODDiff","IntAgeLODFull")]
save_int_age <- arrange(save_int_age, IntAgeChr, IntAgePos)
# save annotated list for sharing
#write.csv(save_int_age, "./QTLscan/output/Threshold6_eQTL_intAge.csv", row.names = FALSE, quote = FALSE)


# Convert transcript and qtl position relative to chromosome positions
# Convert to megabases
chrlen <- sapply(split(Int_age$end, Int_age$chr), max) * 1e-6
chrsum <- c(0, cumsum(chrlen))
names(chrsum) = names(chrlen)

t.gmb <- Int_age$start * 1e-6 # Transcript
q.gmb <- Int_age$IntAgePos * 1e-6 # qtl

# Cumulative sum of previosu positions
for(i in c(2:19, "X")) {

  wh <- which(Int_age$chr == i)
  t.gmb[wh] <- t.gmb[wh] + chrsum[i]

  wh <- which(Int_age$IntAgeChr == i)
    q.gmb[wh] <- q.gmb[wh] + chrsum[i]
}
Int_age$t_gbm <- t.gmb
Int_age$q_gbm <- q.gmb

# Custom lablels & lines
# Only display chr1:19,X
chrtick <- chrsum[1:20]
# Shift axis tick to half way point
max <- max(Int_age$q_gbm)
chrtick_half <- NULL
for (i in 1:length(chrtick)){
  if (i == 20){

    x <- (chrtick[i] + max)/2
    chrtick_half <- c(chrtick_half, x)
  } else {

    x <- (chrtick[i] + chrtick[i + 1])/2
    chrtick_half <- c(chrtick_half, x)
  }
}

# eQTL plot

eQTL <- ggplot(Int_age, aes(x= q_gbm, y= t_gbm)) +
      geom_point(alpha = 0.2) +
      scale_x_continuous("QTL position",
                         breaks = chrtick_half,
                         limits = c(min(Int_age$q_gbm), max(Int_age$q_gbm)),
                         expand = c(0,0)) +
      scale_y_continuous("Gene position",
                         breaks = chrtick_half,
                         limits = c(min(Int_age$t_gbm), max(Int_age$t_gbm)),
                         expand = c(0,0)) +
      geom_vline(xintercept = chrtick[2:20], colour = "grey", size = 0.2) +
      geom_hline(yintercept = chrtick[2:20], colour = "grey", size = 0.2) +
      labs( title = "Interactive-Age eQTLs") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", size = 0.2, fill = NA))

density <- ggplot(Int_age, aes(q_gbm, ..count.., colour = "grey", fill = "grey")) +
      geom_density() +
      scale_colour_manual(name = NA, values = c(grey = "grey"), guide = FALSE) +
      scale_fill_manual(name = NA, values = c(grey = "grey"), guide = FALSE) +
      scale_x_continuous("QTL position",
                         breaks = chrtick_half,
                         limits = c(min(Int_age$q_gbm), max(Int_age$q_gbm)),
                         expand = c(0,0)) +
      scale_y_continuous(name ="Density", breaks = c(0.2, 0.4, 0.6, 0.8)) +
      geom_vline(xintercept = chrtick[2:20], colour = "grey", size = 0.2) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", size = 0.2, fill = NA))

#plot help: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
pdf("./QTLscan/output/plots/eQTL_IntAge_thr8.pdf", width = 9, heigh =9)
eQTL
dev.off()

pdf("./QTLscan/output/plots/eQTL_IntAge_thr8_density.pdf", width = 9, heigh =10)
pushViewport(viewport( layout = grid.layout(10,10)))
print(eQTL, vp = viewport(layout.pos.row = 1:8, layout.pos.col = 1:10))
print(density, vp = viewport(layout.pos.row = 9:10, layout.pos.col = 1:10))
dev.off()
