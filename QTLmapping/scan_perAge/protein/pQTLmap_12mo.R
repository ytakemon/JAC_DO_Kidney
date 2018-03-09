# R/3.4.1
library(tidyverse)
library(grid)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
pQTL_best <- read_csv("./QTLscan/output/12mo_pQTLBestperGene.csv")

# Set LOD threshold 6,8,10
for(thr in c(6,8,10,15)){

print(paste("Preparing threshold @",thr))
LODthreshold <- thr


# Using pQTL_best to create plot
# need to reorder "chr" factors
chr_full <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X","Y", "MT")
AddQTL <- pQTL_best %>%
            mutate(
              chr = factor(chr, levels = chr_full),
              AdditiveChr = factor(AdditiveChr, levels= chr_full)
            ) %>%
            filter(!(chr %in% c("Y","MT"))) %>%
            filter(AdditiveLOD > LODthreshold)

# Annotate postion with genes and save file for sharing
save <- arrange(AddQTL, AdditiveChr, AdditivePos)
# save annotated list for sharing
#write.csv(save, "./QTLscan/output/Threshold6_pQTL_6mo.csv", row.names = FALSE, quote = FALSE)


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

#to adjust pQTL plot a bit to the right to match density
chrtick_halfy <- chrtick_half
names(chrtick_halfy)[20] <- "X  "

# pQTL plot
pQTL <- ggplot(AddQTL, aes(x= q_gbm, y= t_gbm, colour = AdditiveLOD)) +
      geom_point(alpha = 0.8) +
      scale_x_continuous("QTL position",
                         breaks = chrtick_half,
                         limits = c(min(AddQTL$q_gbm), max(AddQTL$q_gbm)),
                         expand = c(0,0)) +
      scale_y_continuous("Gene position",
                         breaks = chrtick_halfy,
                         limits = c(min(AddQTL$t_gbm), max(AddQTL$t_gbm)),
                         expand = c(0,0),
                        sec.axis = dup_axis()) +
      geom_vline(xintercept = chrtick[2:20], colour = "grey", size = 0.2) +
      geom_hline(yintercept = chrtick[2:20], colour = "grey", size = 0.2) +
      labs( title = paste0("12 month pQTLs, threshold @",thr)) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            legend.position = "top",
            panel.border = element_rect(colour = "black", size = 0.2, fill = NA)) +
      scale_colour_gradient2(low = "blue", high = "red", mid = "blue", midpoint = mean(AddQTL$AdditiveLOD))

interval <- seq(0,max(AddQTL$q_gbm), by = 10)
interval <- interval + 5
avgLOD <- NULL
for ( i in interval){
  df <- AddQTL %>% filter((q_gbm > (i-5)) & (q_gbm < (i+5)))
  if(nrow(df) == 0){
    avgLOD <- c(avgLOD,0)
  } else{
    avgLOD <- c(avgLOD, mean(df$AdditiveLOD))
  }
}
df <- data.frame(interval = interval, avgLOD = avgLOD)
z <- c(0,0)
df <- rbind(z, df)

density <- ggplot(AddQTL, aes(q_gbm, colour = "grey", fill = "grey")) +
      geom_histogram(breaks = seq(0,max(AddQTL$q_gbm), by = 10)) +
      geom_line(data = df, aes(x = interval, y = avgLOD*3), colour = "black")+
      scale_colour_manual(name = NA, values = c(grey = "grey"), guide = FALSE) +
      scale_fill_manual(name = NA, values = c(grey = "grey"), guide = FALSE) +
      scale_x_continuous("QTL position",
                         breaks = chrtick_half,
                         limits = c(min(AddQTL$q_gbm), max(AddQTL$q_gbm)),
                         expand = c(0,0)) +
      scale_y_continuous(name ="Density",
                        breaks = seq(0,300, by = 20),
                        sec.axis = sec_axis(trans = ~. / 3, "Agerage LOD")) +
      geom_vline(xintercept = chrtick[2:20], colour = "grey", size = 0.2) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", size = 0.2, fill = NA))

#plot help: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
#pdf("./QTLscan/output/plots/pQTL_Additive_thr8.pdf", width = 9, heigh =9)
#pQTL
#dev.off()
print(paste("Plotting thr @", thr))

pdf(paste0("./QTLscan/output/plots/12mo_pQTL_thr",thr,"_density.pdf"), width = 9, heigh =10)
pushViewport(viewport( layout = grid.layout(10,10)))
print(pQTL, vp = viewport(layout.pos.row = 1:8, layout.pos.col = 1:10))
print(density, vp = viewport(layout.pos.row = 9:10, layout.pos.col = 1:10))
dev.off()
}
