library(qtl)
library(qtlcharts)
source("/hpcdata/cgd/shinyapps/kidney/miniDOQTL.R")
load("/hpcdata/cgd/shinyapps/kidney/shiny_annotation.RData")
load("/hpcdata/cgd/QTL_mapping/kidney_combined/DO192_kidney.RData")

plotType <- 9
gene.name <- "Mrpl20"
level <- "protein"
chr <- 15
file.location <- "/hpcdata/cgd/QTL_mapping/kidney_combined/scanint_protein/ENSMUSP00000030942_Mrpl20.rds"
if (file.exists(file.location)){
  fit <- readRDS(file.location)
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

title <- "Mrpl20"
if (level=="protein") title <- toupper(title)

library(ggplot2)
fit.chr <- fit[fit$chr==chr,]
max.marker <- rownames(fit.chr)[which.max(fit.chr$lod2)[1]]
if (level=="mRNA") ens <- other.ids(gene.name, level)[[1]] else ens <- other.ids(gene.name, level)[[3]]

# if not yet loaded, load the data
if (!exists("expr.mrna")) load("/hpcdata/cgd/QTL_mapping/kidney_combined/DO192_kidney.RData")

if (level=="mRNA") y <- expr.mrna[,ens] else y <- expr.protein[,ens]
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
            ggtitle(max.marker) + mytheme
 print(p)


# get max value of absolute value
dt_max <- dt[which(abs(dt$beta) == max(abs(dt$beta))),]
