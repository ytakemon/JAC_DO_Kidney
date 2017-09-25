library(ggplot2)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
eQTL_add <- read.csv("./QTLscan/output/eQTLAllAdditive.csv")
eQTL_intAge <- read.csv("./QTLscan/output/eQTLAllInteractiveAge.csv")
eQTL_intSex <- read.csv("./QTLscan/output/eQTLAllInteractiveSex.csv")
eQTL_best <- read.csv("./QTLscan/output/eQTLBestperGene.csv")

# Usinb eQTL_best to create plot
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
LODtheshold_additive <- 7.5
LODthreshold_total <- 10.5
LODthreshold_diff <- 7.5

# Plot Interactive-Age eQTLs
# Subset Int-Age LOD above total/full and diff
Int_age <- eQTL_best[(eQTL_best$IntAgeLODDiff > LODthreshold_diff),] # above diff threshold
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

ggplot(Int_age, aes(x= q_gbm, y= t_gbm)) +
      geom_point(alpha = 0.2) +
      scale_x_continuous("QTL position", breaks = ) +
      scale_y_continuous("Gene position") +
      labs( title = "Interactive-Age eQTLs") +
      theme_bw()
