library(dplyr)
library(ggplot2)
library(ggsci)
#library(biomaRt)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO1045_kidney.Rdata")

# Define the following variables
type <- "age" #c("age","add","int")
pheno <- "Alb" #c("Alb","Phs")

# coarse into one column per phenotype in new object: pheno
# DO-0021 is a bug in the cross-sectional that needs to be removed
Upheno_sub <- Upheno %>%
          filter(study == "Cross-sectional" & Mouse.ID != "DO-0021") %>%
          mutate(
            cr = coalesce(cr.u.6, cr.u.12, cr.u.18),
            alb = coalesce(ma.u.6, ma.u.12, ma.u.18),
            phs = coalesce(phs.u.6, phs.u.12, phs.u.18)) %>%
          select(Mouse.ID, cr, alb, phs)

Upheno_sub <- Upheno_sub %>%
              mutate(
                acr = log1p(alb/cr),
                pcr = phs/cr) %>%
              select(Mouse.ID, acr, pcr)

# Set parameters
if (pheno == "Alb"){
  file <- paste0("./QTLscan/addscan_urine/Intscan_alb_all.rds")
  file2 <- paste0("./QTLscan/addscan_urine/Addscan_alb_all.rds")
  chr_select <- "12"
  title <- "log1p(Alb/Cr)"
  ylab <- "log1p(Alb/Cr) (beta +/- SE)"
} else if (pheno == "Phs"){
  file <- paste0("./QTLscan/addscan_urine/Intscan_phs_all.rds")
  file2 <- paste0("./QTLscan/addscan_urine/Addscan_phs_all.rds")
  chr_select <- "2"
  title <- "Phs/Cr"
  ylab <- "Phs/Cr (beta +/- SE)"
}

# Get LOD scores
if(type == "age"){
  intfit <- readRDS(file)
  addfit <- readRDS(file2)
  fit <- intfit - addfit
} else if(type == "int"){
  fit <- readRDS(file)
} else if(type == "add"){
  fit <- readRDS(fitl2)
}
fit <- as.data.frame(fit)

# Max Marker
fit$chr <- MM_snps[MM_snps$marker %in% rownames(fit),]$chr
fit.chr <- fit[fit$chr==chr_select,]
max.marker <- rownames(fit.chr)[which.max(fit.chr[,1])]
max.pos <- MM_snps[MM_snps$marker == max.marker,]$pos

# Subset data
Upheno_sub <- bind_cols(Upheno_sub,
              samples %>%
              filter(Mouse.ID %in% Upheno_sub$Mouse.ID) %>%
              select(Sex, Cohort.Age.mo)) %>%
              filter(!is.na(Cohort.Age.mo))
# Subset genoprobs
probs <- genoprobs[Upheno_sub$Mouse.ID,,]
# calculate divation form mean
if(pheno == "Alb") y <- Upheno_sub$acr else y <- Upheno_sub$pcr
y[Upheno_sub$Sex == "M"] <- y[Upheno_sub$Sex == "M"] - mean(y[Upheno_sub$Sex == "M"], na.rm = TRUE)
y[Upheno_sub$Sex == "F"] <- y[Upheno_sub$Sex == "F"] - mean(y[Upheno_sub$Sex == "F"], na.rm = TRUE)

# Get effect by age:
coef <- se <- tstat <- NULL
for( age in c("6", "12", "18")){
  sel <- Upheno_sub$Cohort.Age.mo == age
  lm <- summary(lm(y[sel] ~ 0 + probs[sel,,max.marker]))
  coef <- bind_cols(coef, as.data.frame(lm$coef[1:8,1]))
  se <- bind_cols(se, as.data.frame(lm$coef[1:8,2]))
  tstat <- bind_cols(tstat, as.data.frame(lm$coef[1:8,3]))
}

# Compile into df
dt <- data.frame(Allele = LETTERS[rep(1:8,3)],
           Age = factor(rep(c(6,12,18), each=8)),
           beta = unlist(coef),
           betase = unlist(se),
           Tstat = unlist(tstat))

pdf(paste0("./QTLscan/output/plots/Urine_",pheno,"_AlleleEffect_byAge_chr", chr_select,".pdf"), width = 12, height = 6)
ggplot(dt, aes(x = as.numeric(Allele), y = beta, colour = Age)) +
      geom_errorbar(aes(ymin=beta - betase, ymax=beta + betase), width=.3, position = position_dodge(0.3)) +
      geom_point(aes(size=abs(Tstat)), position = position_dodge(0.3)) +
      xlab("Allele") +
      ylab( ylab ) +
      scale_x_continuous(breaks=c(1:8), labels=LETTERS[1:8]) +
      geom_abline(intercept = 0, slope = 0, colour=I("grey")) +
      ggtitle(paste0(title, " : ", max.marker, " (Chr",chr_select,":",max.pos, " Mbp)")) +
      theme(panel.border = element_rect(linetype = "solid", fill=NA, colour = "grey")) +
      scale_color_aaas()
dev.off()
