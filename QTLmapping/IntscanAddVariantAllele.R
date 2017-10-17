library(qtl)
library(qtlcharts)
library(ggplot2)
source("/hpcdata/cgd/shinyapps/kidney/miniDOQTL.R")
load("/hpcdata/cgd/shinyapps/kidney/shiny_annotation.RData")
load("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/RNAseq_data/DO188b_kidney.RData")
eQTL <- read.csv("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/QTLscan/output/Threshold8_eQTL_intAge.csv")
pQTL <- read.csv("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/QTLscan/output/Threshold8_pQTL_intAge.csv")

# eQTL variant----------------------------------------------------------------
eQTL$AlleleVariant <- NA
eQTL$AgeVariant <- NA
for (i in 1:nrow(eQTL)){
  gene.name <- eQTL$id[i]
  symbol <- eQTL$symbol[i]
  level <- "mRNA"
  chr <- eQTL$chr[i]
  file <- paste0("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/QTLscan/intscan_mrna/Age/",
                 gene.name,
                 "_",
                 symbol,
                 ".rds")
  # list.files(path = "/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/QTLscan/intscan_mrna/Age", pattern = "gene.name")
  if (file.exists(file)){
    fit <- readRDS(file)
  } else {
    next
  }

  # for given MGI symbol, find Ensembl ids
  other.ids <- function(gene.name, level) {
    if (level == "mRNA") {
      sel <- which(mRNA.list$symbol == gene.name)[1]
      if (!is.na(sel)) return(mRNA.list[sel,]) else return(c(NA,NA,NA))
    }
    if (level == "protein") {
      sel <- which(protein.list$symbol == gene.name)[1]
      if (!is.na(sel)) return(protein.list[sel,]) else return(c(NA,NA,NA))
    }
  }

  fit <- as.data.frame(fit)
  fit$chr <- substr(rownames(fit),1,regexpr("_",rownames(fit)) - 1)
  fit.chr <- fit[fit$chr==chr,]
  max.marker <- rownames(fit.chr)[which.max(fit.chr$pheno1)]
  if (level=="mRNA") ens <- other.ids(symbol, level)[[1]] else ens <- other.ids(symbol, level)[[3]]
  if (level=="mRNA") y <- expr.mrna[,ens] else y <- expr.protein[,ens]
  y[annot.samples$Sex=="M"] <- y[annot.samples$Sex=="M"] - mean(y[annot.samples$Sex=="M"], na.rm=TRUE)
  y[annot.samples$Sex=="F"] <- y[annot.samples$Sex=="F"] - mean(y[annot.samples$Sex=="F"], na.rm=TRUE)
  coef <- se <- tstat <-  NULL

  for (age in c("6", "12", "18")) {
    sel <- annot.samples$Age == age
    lm.fit <- summary(lm(y[sel] ~ 0 + genoprobs[sel,,max.marker]))
    coef <- cbind(coef, lm.fit$coef[1:8,1])
    se   <- cbind(se, lm.fit$coef[1:8,2])
    tstat   <- cbind(tstat, lm.fit$coef[1:8,3])
  }

  dt <- data.frame(Allele = LETTERS[rep(1:8,3)],
             Age = factor(rep(c(6,12,18), each=8)),
             beta = as.vector(coef),
             beta.se = as.vector(se),
             Tstat=as.vector(tstat))

  dt_max <- dt[which(abs(dt$beta) == max(abs(dt$beta))),]
  eQTL$AlleleVariant[i] <- as.character(dt_max$Allele[1])
  eQTL$AgeVariant[i] <- as.character(dt_max$Age[1])
}
write.csv(eQTL, file = "/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/QTLscan/output/Threshold8_eQTL_intAge_annotated.csv")

# Plotting -------------------------------------------------------------------
# mytheme <- theme(panel.border = element_rect(linetype = "solid", , fill=NA, colour = "grey"))
# p <- ggplot(aes(x=as.numeric(Allele)+as.numeric(Age)/15-2/15,
#                 y=beta,colour=Age), data=dt) +
#             geom_errorbar(aes(ymin=beta-beta.se, ymax=beta+beta.se), width=.1) +
#             geom_point(aes(size=abs(Tstat))) + xlab("Allele") +
#             ylab("beta +/- SE") +
#             scale_x_continuous(breaks=c(1:8), labels=LETTERS[1:8]) +
#             geom_abline(intercept = 0, slope = 0, colour=I("grey")) +
#             ggtitle(max.marker) + mytheme
# print(p)

# get max value of absolute value ---------------------------------------------
