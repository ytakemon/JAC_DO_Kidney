# Is age/sex mediated by RNA? (p-value version)
# See Petr's page: http://simecek.xyz/TheAgingProteome/new-med-plots.html

# References:
# https://github.com/simecek/TheAgingProteome/blob/master/analysis/new-med-plots.Rmd

# Usage on Cadillac
# cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts
# qsub -v I=kidney_anova_fdr_output.csv,script=Mediation_analysis Rsubmit_args.sh

# R/3.3.2
library(ggplot2)
library(dplyr)
library(scales)

# Load data
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
args <- commandArgs(trailingOnly = TRUE)
file <- list.files(pattern = paste0("^",args[[1]]), recursive = TRUE)

# -log10(p) transformation for ggplot
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
              log_breaks(base = base),
              domain = c(1e-100, Inf))
}

# read ANOVA table - kidney
dt <- read.csv(file[[1]], as.is = TRUE) %>% select(symbol, starts_with("p.Prot_Age"),
         starts_with("p.Prot_Sex"))

# gathering `dt` from 4 cols to 2 cols (->ggplot)
# similar to the melt() process from reshape2 package.
tmp1 <- select(dt, symbol, starts_with("p.Prot_Age"))
tmp2 <- select(dt, symbol, starts_with("p.Prot_Sex"))
names(tmp1) <- names(tmp2) <- c("symbol", "x", "y")
tmp1$var <- "Age"
tmp2$var <- "Sex"
dt2 <- rbind(tmp1, tmp2)

pl <- ggplot(dt2, aes(x=x,  y=y, text=symbol)) +
  geom_point(alpha=0.2) +
  geom_abline(intercept = 0, slope = 1, color="red") +
  scale_x_continuous(trans=reverselog_trans(10)) +
  scale_y_continuous(trans=reverselog_trans(10)) +
  facet_wrap(~var) +
  xlab("p-value (X)") +
  ylab("p-value (X | mRNA)") +
  theme_bw() +
  labs(title="Kidney")

pdf("./Plot/mRNA_mediation.pdf", width = 12, height = 6)
pl
dev.off()
