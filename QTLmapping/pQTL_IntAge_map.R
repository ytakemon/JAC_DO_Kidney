# R/3.4.1
library(ggplot2)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
pQTL_add <- read.csv("./QTLscan/output/pQTLAllAdditive.csv")
pQTL_intAge <- read.csv("./QTLscan/output/pQTLAllInteractiveAge.csv")
pQTL_intSex <- read.csv("./QTLscan/output/pQTLAllInteractiveSex.csv")
pQTL_best <- read.csv("./QTLscan/output/pQTLBestperGene.csv")

# Usinb pQTL_best to create plot
# need to reorder "chr" factors
chr_full <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X","Y", "MT")
pQTL_best$chr <- factor(pQTL_best$chr, levels = chr_full)
pQTL_best$IntAgeChr <- factor(pQTL_best$IntAgeChr, levels= chr_full)
pQTL_best$IntSexChr <- factor(pQTL_best$IntSexChr, levels= chr_full)
# Subset out chr 1-19,X from data
chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X")
pQTL_best <- pQTL_best[pQTL_best$chr %in% chr, ]

# Combine chr&positions into one numerical verctor
pQTL_best$Gene_point <- paste(pQTL_best$chr, pQTL_best$start, sep = ".")
pQTL_best$IntAge_point <- paste(pQTL_best$IntAgeChr, pQTL_best$IntAgePos, sep = ".")
pQTL_best$IngSex_point <- paste(pQTL_best$IntSexChr, pQTL_best$IntSexPos, sep = ".")

# Set LOD threshold
LODtheshold_additive <- 7.5
LODthreshold_total <- 10.5
LODthreshold_diff <- 7.5

# Plot Interactive-Age pQTLs
# Subset Int-Age LOD above total/full and diff
Int_age <- pQTL_best[(pQTL_best$IntAgeLODDiff > LODthreshold_diff),] # above diff threshold
Int_age <- Int_age[(Int_age$IntAgeLODFull > LODthreshold_total),] # above full/total threshold

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
# Shift thick to half way point
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

# pQTL plot
pdf("./QTLscan/output/plots/pQTL_IntAge.pdf", width = 6, heigh =6)
ggplot(Int_age, aes(x= q_gbm, y= t_gbm)) +
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
      labs( title = "Interactive-Age pQTLs") +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", size = 0.2, fill = NA))
dev.off()
