# Leaving this for discussion, incomplete
library(dplyr)
library(ggsci)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney.RData")

Markers <- annot.samples
Markers$Age <- as.factor(as.character(Markers$Age))
Markers$Age <- factor(Markers$Age, levels(Markers$Age)[c(3,1,2)]) # Change levels order

# find gene id
Havcr1_id <- annot.mrna[annot.mrna$symbol == "Havcr1",]$id
Lcn2_id <- annot.mrna[annot.mrna$symbol == "Lcn2",]$id
Naglu_id <- annot.mrna[annot.mrna$symbol == "Naglu",]$id
Cst3_id <- annot.mrna[annot.mrna$symbol == "Cst3",]$id

#Havcr1_pid <- annot.protein[annot.protein$symbol == "Havcr1",]$id doesn't exist
Lcn2_pid <- annot.protein[annot.protein$symbol == "Lcn2",]$id
Naglu_pid <- annot.protein[annot.protein$symbol == "Naglu",]$id
Cst3_pid <- annot.protein[annot.protein$symbol == "Cst3",]$id

Markers$Havcr1_rna_raw <- raw.mrna[,Havcr1_id]
#Markers$Havcr1_prot_raw <- raw.protein[,Havcr1_pid] doesn't exist
Markers$Lcn2_rna_raw <- raw.mrna[,Lcn2_id]
Markers$Lcn2_prot_raw <- raw.protein[,Lcn2_pid]
Markers$Naglu_rna_raw <- raw.mrna[,Naglu_id]
Markers$Naglu_prot_raw <- raw.protein[,Naglu_pid]
Markers$Cst3_rna_raw <- raw.mrna[,Cst3_id]
Markers$Cst3_prot_raw <- raw.protein[,Cst3_pid]



#pdf("./Plot/Akt1Allele_Akt1exp_by_Age.pdf", width = 12, height = 7)
ggplot(Markers, aes(x = Age, y = Havcr1_rna_raw, colour = Age)) +
  geom_smooth(method = "lm", se = FALSE, aes(group = Age, colour = Age)) +
  theme_bw() +
  labs( title = "Akt1 mRNA expression by genotype",
        subtitle = paste0("Female: 93, (NZO/NZO=0 , NZO/Other=9 , Other/Other=84 ) \nMale: 95, (NZO/NZO=3 , NZO/Other=12 , Other/Other=80 )"),
        y = "Akt1 mRNA expression",
        x = "Age") +
  facet_grid(. ~ Sex) +
  scale_colour_aaas()
#dev.off()
