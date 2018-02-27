library(dplyr)
library(ggplot2)
library(ggsci)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")
load("./shiny_annotation.RData")

gene_select <- "Qrsl1"
level_select <- "protein"
chr <- 7

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

query <- other.ids(gene_select, level_select)
# Gather gene info for query
if (level_select == "mRNA"){
  id <- query$id[1]
  symbol <- query$symbol[1]
  file <- paste0("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/QTLscan/intscan_mrna/Age/",
                 id,
                 "_",
                 symbol,
                 ".rds")
} else {
  id <- query$protein_id[1]
  symbol <- query$symbol[1]
  file <- paste0("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/QTLscan/intscan_prot_pbatch/Age/",
                 id,
                 "_",
                 symbol,
                 ".rds")
}

# Validate file existance
if (file.exists(file)){
  fit <- readRDS(file)
} else {
  fit <- NA
}

# fit
fit <- as.data.frame(fit)
fit$chr <- substr(rownames(fit),1,regexpr("_",rownames(fit)) - 1)
fit.chr <- fit[fit$chr==chr,]
max.marker <- rownames(fit.chr)[which.max(fit.chr[,1])[1]]
max.marker.pos <- snps[snps$marker == max.marker, ]$bp

ens <- other.ids(gene.name, level)[[1]]
if (level_select=="mRNA") ens <- other.ids(symbol, level_select)[[1]] else ens <- other.ids(symbol, level_select)[[3]]

# Get data
if (level_select=="mRNA") y <- expr.mrna[,ens] else y <- expr.protein[,ens]
y[annot.samples$Sex=="M"] <- y[annot.samples$Sex=="M"] - mean(y[annot.samples$Sex=="M"], na.rm=TRUE)
y[annot.samples$Sex=="F"] <- y[annot.samples$Sex=="F"] - mean(y[annot.samples$Sex=="F"], na.rm=TRUE)

# Get effect by age:
coef <- se <- tstat <- NULL
for( age in c("6", "12", "18")){
  sel <- annot.samples$Age == age
  lm <- summary(lm(y[sel] ~ 0 + genoprobs[sel,,max.marker]))
  coef <- cbind(coef, lm$coef[1:8,1])
  se   <- cbind(se, lm$coef[1:8,2])
  tstat   <- cbind(tstat, lm$coef[1:8,3])
}

# Compile into df
dt <- data.frame(Allele = LETTERS[rep(1:8,3)],
           Age = factor(rep(c(6,12,18), each=8)),
           beta = as.vector(coef),
           betase = as.vector(se),
           Tstat=as.vector(tstat))

pdf(paste0("./QTLscan/output/plots/AlleleEffect_intQTL_",level_select,"_",gene_select,".pdf"),  width = 12, height = 6)
ggplot(dt, aes(x = as.numeric(Allele), y = beta, colour = Age)) +
  geom_errorbar(aes(ymin=beta - betase, ymax=beta + betase), width=.3, position = position_dodge(0.3)) +
  geom_point(aes(size=abs(Tstat)), position = position_dodge(0.3)) +
  ylab("beta +/- SE") +
  xlab("Alleles") +
  scale_x_continuous(breaks=c(1:8), labels=LETTERS[1:8]) +
  geom_abline(intercept = 0, slope = 0, colour=I("grey")) +
  ggtitle(paste0(gene_select, " intQTL \nMax marker",max.marker, " (Chr",chr,":",max.marker.pos, " Mbp)")) +
  theme(panel.border = element_rect(linetype = "solid", fill=NA, colour = "grey")) +
  scale_color_aaas()
dev.off()
