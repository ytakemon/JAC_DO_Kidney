# R/3.4.1
library(tidyverse)
library(grid)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
pQTL_best <- read.csv("./QTLscan/output/pQTLAllAdditive_pbatch.csv")

# Using pQTL_best to create plot
# need to reorder "chr" factors
chr_full <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X","Y", "MT")
pQTL_best$chr <- factor(pQTL_best$chr, levels = chr_full)
pQTL_best$AdditiveChr <- factor(pQTL_best$AdditiveChr, levels= chr_full)

# Subset out chr 1-19,X from data
chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X")
pQTL_best <- pQTL_best[pQTL_best$chr %in% chr, ]

# Combine chr&positions into one numerical verctor
pQTL_best$Gene_point <- paste(pQTL_best$chr, pQTL_best$start, sep = ".")
pQTL_best$Add_point <- paste(pQTL_best$AdditiveChr, pQTL_best$AdditivePos, sep = ".")


# Set LOD threshold
LODthreshold <- 8.5

# Plot Interactive-Age pQTLs
# Subset Int-Age LOD above total/full and diff
AddQTL <- pQTL_best[(pQTL_best$AdditiveLOD > LODthreshold),] # above diff threshold

# Convert transcript and qtl position relative to chromosome positions
# Convert to megabases
chrlen <- sapply(split(AddQTL$end, AddQTL$chr), max) * 1e-6
chrsum <- c(0, cumsum(chrlen))
names(chrsum) = names(chrlen)

t.gmb <- AddQTL$start * 1e-6 # Transcript
q.gmb <- AddQTL$AdditivePos * 1e-6 # qtl

# Cumulative sum of previosu positions
for(i in c(2:19, "X")) {

  wh <- which(AddQTL$chr == i)
  t.gmb[wh] <- t.gmb[wh] + chrsum[i]

  wh <- which(AddQTL$AdditiveChr == i)
    q.gmb[wh] <- q.gmb[wh] + chrsum[i]
}
AddQTL$t_gbm <- t.gmb
AddQTL$q_gbm <- q.gmb

# Custom lablels & lines
# Only display chr1:19,X
chrtick <- chrsum[1:20]
# Shift axis tick to half way point
max <- max(AddQTL$q_gbm)
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

mPlot <-  ggplot(AddQTL, aes(x= q_gbm, y= t_gbm)) +
      geom_point(alpha = 0.2) +
      scale_x_continuous("QTL position",
                         breaks = chrtick_half,
                         limits = c(min(AddQTL$q_gbm), max(AddQTL$q_gbm)),
                         expand = c(0,0)) +
      scale_y_continuous("Gene position",
                         breaks = chrtick_half,
                         limits = c(min(AddQTL$t_gbm), max(AddQTL$t_gbm)),
                         expand = c(0,0)) +
      geom_vline(xintercept = chrtick[2:20], colour = "grey", size = 0.2) +
      geom_hline(yintercept = chrtick[2:20], colour = "grey", size = 0.2) +
      labs( title = "Additive pQTLs (LOD threshold > 8.5)") +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", size = 0.2, fill = NA))

density <- ggplot(AddQTL, aes(q_gbm, colour = "grey", fill = "grey")) +
      geom_histogram(breaks = seq(0,max(AddQTL$q_gbm), by = 10)) +
      scale_colour_manual(name = NA, values = c(grey = "grey"), guide = FALSE) +
      scale_fill_manual(name = NA, values = c(grey = "grey"), guide = FALSE) +
      scale_x_continuous("QTL position",
                         breaks = chrtick_half,
                         limits = c(min(AddQTL$q_gbm), max(AddQTL$q_gbm)),
                         expand = c(0,0)) +
      scale_y_continuous(name ="Density",
                         breaks = seq(0,450, by = 20)) +
      geom_vline(xintercept = chrtick[2:20], colour = "grey", size = 0.2) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", size = 0.2, fill = NA))

# pQTL plot
pdf("./QTLscan/output/plots/Add_pQTL_thr8.5.pdf", width = 9, height =10)
pushViewport(viewport( layout = grid.layout(10,10)))
print(mPlot, vp = viewport(layout.pos.row = 1:8, layout.pos.col = 1:10))
print(density, vp = viewport(layout.pos.row = 9:10, layout.pos.col = 1:10))
dev.off()
