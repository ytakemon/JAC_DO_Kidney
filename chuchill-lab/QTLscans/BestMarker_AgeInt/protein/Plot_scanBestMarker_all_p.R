# R/3.4.1
library(tidyverse)
library(grid)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
best <- read.csv("./QTLscan/scanBestMarker_protein/BestMarker_BestperGene_protein.csv")

# check distribution of LOD scores
hist(best$IntAgeLODDiff, breaks =100)
# picking blunt threshold
LODthreshold_diff <- 7.5

# Using best to create plot
# need to reorder "chr" factors
chr_full <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X","Y", "MT")
best$chr <- factor(best$chr, levels = chr_full)
best$IntAgeChr <- factor(best$IntAgeChr, levels= chr_full)

# Subset out chr 1-19,X from data
chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X")
best <- best[best$chr %in% chr, ]

# Plot Interactive-Age eQTLs
# Subset Int-Age LOD above total/full and diff
Int_age <- best[(best$IntAgeLODDiff > LODthreshold_diff),] # above diff threshold

# Annotate Interactive-Age postion with genes and save file for sharing
#save_int_age <- arrange(Int_age, IntAgeChr, IntAgePos)
#write_csv(save_int_age, path = paste0("./QTLscan/scanBestMarker_protein/BestMarker_BestperGene_protein_thr",LODthreshold_diff,".csv"))


# Convert transcript and qtl position relative to chromosome positions
# Convert to megabases
chrlen <- sapply(split(Int_age$end, Int_age$chr), max) * 1e-6
chrsum <- c(0, cumsum(chrlen))
names(chrsum) = names(chrlen)

t.gmb <- Int_age$start * 1e-6 # Transcript
q.gmb <- Int_age$IntAgePos # qtl

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

#to adjust eQTL plot a bit to the right to match density
chrtick_halfy <- chrtick_half
names(chrtick_halfy)[20] <- "X  "

# eQTL plot
mPlot <- ggplot(Int_age, aes(x= q_gbm, y= t_gbm)) +
      geom_point(alpha = 0.2) +
      scale_x_continuous("QTL position",
                         breaks = chrtick_half,
                         limits = c(min(Int_age$q_gbm), max(Int_age$q_gbm)),
                         expand = c(0,0)) +
      scale_y_continuous("Gene position",
                         breaks = chrtick_halfy,
                         limits = c(min(Int_age$t_gbm), max(Int_age$t_gbm)),
                         expand = c(0,0)) +
      geom_vline(xintercept = chrtick[2:20], colour = "grey", size = 0.2) +
      geom_hline(yintercept = chrtick[2:20], colour = "grey", size = 0.2) +
      labs( title = "Interactive-Age pQTLscan by Marker") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            legend.position = "top",
            panel.border = element_rect(colour = "black", size = 0.2, fill = NA))

interval <- seq(0,max(Int_age$q_gbm), by = 10)
interval <- interval + 5
avgLOD <- NULL
for ( i in interval){
  df <- Int_age %>% filter((q_gbm > (i-5)) & (q_gbm < (i+5)))
  if(nrow(df) == 0){
    avgLOD <- c(avgLOD,0)
  } else{
    avgLOD <- c(avgLOD, mean(df$IntAgeLODDiff))
  }
}
df <- data.frame(interval = interval, avgLOD = avgLOD)
z <- c(0,0)
df <- rbind(z, df)

density <- ggplot(Int_age, aes(q_gbm, colour = "grey", fill = "grey")) +
      geom_histogram(breaks = seq(0,max(Int_age$q_gbm), by = 10)) +
      scale_colour_manual(name = NA, values = c(grey = "grey"), guide = FALSE) +
      scale_fill_manual(name = NA, values = c(grey = "grey"), guide = FALSE) +
      scale_x_continuous("QTL position",
                         breaks = chrtick_half,
                         limits = c(min(Int_age$q_gbm), max(Int_age$q_gbm)),
                         expand = c(0,0)) +
      scale_y_continuous(name ="Density",
                         breaks = seq(0,450, by = 10)) +
      geom_vline(xintercept = chrtick[2:20], colour = "grey", size = 0.2) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", size = 0.2, fill = NA))

#plot help: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
pdf(paste0("./QTLscan/scanBestMarker_protein/BestMarker_BestperGene_protein_thr",LODthreshold_diff,".pdf"), width = 9, heigh =10)
pushViewport(viewport( layout = grid.layout(10,10)))
print(mPlot, vp = viewport(layout.pos.row = 1:8, layout.pos.col = 1:10))
print(density, vp = viewport(layout.pos.row = 9:10, layout.pos.col = 1:10))
dev.off()
