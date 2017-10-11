library(shiny)
library(qtl)
library(qtlcharts)
library(ggplot2)
source("/hpcdata/cgd/shinyapps/kidney/miniDOQTL.R")
load("/hpcdata/cgd/shinyapps/kidney/shiny_annotation.RData")
load("/hpcdata/cgd/QTL_mapping/kidney_combined/DO192_kidney.RData")

plotType <- 9
gene.name <- "Slc6a20b"
level <- "mRNA"
chr <- 9

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

m.file.location <- paste0("/hpcdata/cgd/QTL_mapping/kidney_combined/scanint_mrna/",other.ids(gene.name,level)[[1]],"_",gene.name,".rds")
p.file.location <- paste0("/hpcdata/cgd/QTL_mapping/kidney_combined/scanint_protein/",other.ids(gene.name,level)[[3]],"_",gene.name,".rds")

if (file.exists(m.file.location)){
  m.fit <- readRDS(m.file.location)
}
if (file.exists(p.file.location)){
  p.fit <- readRDS(p.file.location)
}

# mRNA ------------------------------------------------------------------------
m.fit.chr <- m.fit[m.fit$chr==chr,]
max.marker <- rownames(m.fit.chr)[which.max(m.fit.chr$lod2)[1]]
ens <- other.ids(gene.name, level)[[1]]
y <- expr.mrna[, ens]
y[pheno$Sex=="M"] <- y[pheno$Sex=="M"] - mean(y[pheno$Sex=="M"], na.rm=TRUE)
y[pheno$Sex=="F"] <- y[pheno$Sex=="F"] - mean(y[pheno$Sex=="F"], na.rm=TRUE)
coef <- se <- tstat <-  NULL
for (age in c("6", "12", "18")) {
  sel <- pheno$Age == age
  lm.fit <- summary(lm(y[sel] ~ 0 + probs[sel,,max.marker]))
  coef <- cbind(coef, lm.fit$coef[1:8,1])
  se   <- cbind(se, lm.fit$coef[1:8,2])
  tstat   <- cbind(tstat, lm.fit$coef[1:8,3])
}
dt <- data.frame(Allele = LETTERS[rep(1:8,3)],
           Age = factor(rep(c(6,12,18), each=8)),
           beta = as.vector(coef),
           beta.se = as.vector(se),
           Tstat=as.vector(tstat))
mytheme <- theme(panel.border = element_rect(linetype = "solid", , fill=NA, colour = "grey"))
p <- ggplot(aes(x=as.numeric(Allele)+as.numeric(Age)/15-2/15,
                y=beta,colour=Age), data=dt) +
            geom_errorbar(aes(ymin=beta-beta.se, ymax=beta+beta.se), width=.1) +
            geom_point(aes(size=abs(Tstat))) + xlab("Allele") +
            ylab("beta +/- SE") +
            scale_x_continuous(breaks=c(1:8), labels=LETTERS[1:8]) +
            geom_abline(intercept = 0, slope = 0, colour=I("grey")) +
            ggtitle(paste0(max.marker, " ", gene.name, " allele effect by age")) + mytheme

pdf(paste0("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Plot/AlleleEffect/AlleleEffectRNA_",gene.name,".pdf"), height = 7, width = 7)
print(p)
dev.off()

# Protein ----------------------------------------------------------------------
title <- gene.name
p.title <- toupper(title)
p.fit.chr <- p.fit[p.fit$chr==chr,]
max.marker <- rownames(p.fit.chr)[which.max(p.fit.chr$lod2)[1]]
ens <- other.ids(gene.name, level)[[3]]
y <- expr.protein[, ens]
y[pheno$Sex=="M"] <- y[pheno$Sex=="M"] - mean(y[pheno$Sex=="M"], na.rm=TRUE)
y[pheno$Sex=="F"] <- y[pheno$Sex=="F"] - mean(y[pheno$Sex=="F"], na.rm=TRUE)
coef <- se <- tstat <-  NULL
for (age in c("6", "12", "18")) {
  sel <- pheno$Age == age
  lm.fit <- summary(lm(y[sel] ~ 0 + probs[sel,,max.marker]))
  coef <- cbind(coef, lm.fit$coef[1:8,1])
  se   <- cbind(se, lm.fit$coef[1:8,2])
  tstat   <- cbind(tstat, lm.fit$coef[1:8,3])
}

dt <- data.frame(Allele = LETTERS[rep(1:8,3)],
           Age = factor(rep(c(6,12,18), each=8)),
           beta = as.vector(coef),
           beta.se = as.vector(se),
           Tstat=as.vector(tstat))
mytheme <- theme(panel.border = element_rect(linetype = "solid", , fill=NA, colour = "grey"))
p <- ggplot(aes(x=as.numeric(Allele)+as.numeric(Age)/15-2/15,
                y=beta,colour=Age), data=dt) +
            geom_errorbar(aes(ymin=beta-beta.se, ymax=beta+beta.se), width=.1) +
            geom_point(aes(size=abs(Tstat))) + xlab("Allele") +
            ylab("beta +/- SE") +
            scale_x_continuous(breaks=c(1:8), labels=LETTERS[1:8]) +
            geom_abline(intercept = 0, slope = 0, colour=I("grey")) +
            ggtitle(paste0(max.marker, " ", p.title, " allele effect by age")) + mytheme

pdf(paste0("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Plot/AlleleEffect/AlleleEffectProt_",p.title,".pdf"), height = 7, width = 7)
print(p)
dev.off()


# get max value of absolute value
dt_max <- dt[which(abs(dt$beta) == max(abs(dt$beta))),]
