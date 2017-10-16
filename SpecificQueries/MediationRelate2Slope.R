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
file <- "./Anova_output/kidney_anova_output.csv"
dt <- read.csv(file[[1]], as.is = TRUE) %>% select(id, symbol, starts_with("p.Prot_Age"),
         starts_with("p.Prot_Sex"))
file <- "./Anova_output/gene_lists/Quad_Age_Sig.csv"
slope <- read.csv(file[[1]], as.is = TRUE) %>% select(id,symbol, starts_with("quad"))
dt_slope <- dt[dt$id %in% slope$id,]
dt_slope$quadI <- slope$quadI
dt_slope$quadII <- slope$quadII
dt_slope$quadIII <- slope$quadIII
dt_slope$quadIV <- slope$quadIV

quadI <- dt_slope[dt_slope$quadI == TRUE, ]
quadII <- dt_slope[dt_slope$quadII == TRUE, ]
quadIII <- dt_slope[dt_slope$quadIII == TRUE, ]
quadIV <- dt_slope[dt_slope$quadIV == TRUE, ]

# -log10(p) transformation for ggplot
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
              log_breaks(base = base),
              domain = c(1e-100, Inf))
}

# read ANOVA table - kidney
#


# gathering `dt` from 4 cols to 2 cols (->ggplot)
# similar to the melt() process from reshape2 package.
Age <- select(dt, symbol, starts_with("p.Prot_Age"))
names(Age) <- c("symbol", "x", "y")
Age$var <- "Age"

pdf("./Plot/MediationRelate2Slope.pdf", width = 12, height = 6)
ggplot() +
  geom_point(data = Age, aes(x = x, y = y, colour = "total"), alpha=0.5) +
  geom_point(data = quadI, aes(x = p.Prot_Age.Sex, y = p.Prot_Age.SexmRNA, colour = "quadI")) +
  geom_point(data = quadII, aes(x = p.Prot_Age.Sex, y = p.Prot_Age.SexmRNA, colour = "quadII")) +
  geom_point(data = quadIII, aes(x = p.Prot_Age.Sex, y = p.Prot_Age.SexmRNA, colour = "quadIII")) +
  geom_point(data = quadIV, aes(x = p.Prot_Age.Sex, y = p.Prot_Age.SexmRNA, colour = "quadIV")) +
  geom_abline(intercept = 0, slope = 1, color="red") +
  scale_colour_manual(name = "Quadrants", values = c(total = "grey", quadI = "red", quadII = "blue", quadIII = "green", quadIV = "orange")) +
  scale_x_continuous(trans=reverselog_trans(10)) +
  scale_y_continuous(trans=reverselog_trans(10)) +
  xlab("p-value (X)") +
  ylab("p-value (X | mRNA)") +
  theme_bw() +
  labs(title="Kidney ~ Age")
dev.off()
