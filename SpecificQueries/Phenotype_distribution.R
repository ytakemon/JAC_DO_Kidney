library(dplyr)
library(reshape)
library(ggplot2)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO1045_kidney.Rdata")

# Subset pheontype: Alb
pheno <- Upheno[Upheno$study == "Cross-sectional",]
pheno <- pheno[-1,]
pheno_cr <- pheno[,c("Mouse.ID", "cr.u.6", "cr.u.12", "cr.u.18")]
pheno_ma <- pheno[,c("Mouse.ID", "ma.u.6", "ma.u.12", "ma.u.18")]

# colapse to one column and combine
pheno_cr <- pheno_cr %>% mutate(
              cr.u.all = coalesce(cr.u.6, cr.u.12, cr.u.18)
            )
pheno_ma <- pheno_ma %>% mutate(
              ma.u.all = coalesce(ma.u.6, ma.u.12, ma.u.18)
            )
pheno <- pheno_cr[,c("Mouse.ID","cr.u.all")]
pheno$ma.u.all <- pheno_ma$ma.u.all
pheno <- pheno[!is.na(pheno[,2] & pheno[,3]),]
pheno$ACR <- pheno$ma.u.all / pheno$cr.u.all
rownames(pheno) <- pheno$Mouse.ID

# Subset dataset
samples <- samples[pheno$Mouse.ID,]

# Get age
pheno$age <- samples$Cohort.Age.mo
pheno$sex <- samples$Sex
pheno$age <- as.factor(as.character(pheno$age))
phenoF <- pheno[pheno$sex == "F",]


ggplot(phenoF, aes(x = log(ACR), y =..count.., fill = age, colour = age)) +
  geom_density( alpha = 0.1)+
  scale_x_continuous("Albumin to creatinine ratio (mg/g)") +
  scale_y_continuous(" ") +
  labs( title = "Females") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  coord_flip()
pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4a5_fig2_ACR_GFR_dist_YUKA.pdf", width = 12, height = 6)
pushViewport(viewport(layout = grid.layout(1, 2)))
print(ggplot_F, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(ggplot_M, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()
