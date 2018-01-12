# Identify animals that have the average given expression values for the given group.
library(dplyr)
library(ggplot2)
library(ggsci)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")
load("./shiny_annotation.RData")

# Selection quality:
Gene_name <- "Lrp2"
Data <- "protein"
Sex <- "F"
Age <- c(6, 18)
Num <- 5

samples <- annot.samples

# Identify gene name
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

gene <- other.ids(Gene_name, Data)
# Get expression data
samples$RNAexpr <- expr.mrna[, gene$id]
samples$Protexpr <- expr.protein[, gene$protein_id]

# Subset selction
samples <- samples[ samples$Age %in% Age,]
samples <- samples[ samples$Sex %in% Sex,]
cal <- aggregate(Protexpr ~ Age, data = samples, mean)
mo6 <- samples[ samples$Age == 6,] %>% arrange(Protexpr)
mo18 <- samples[ samples$Age == 18,] %>% arrange(Protexpr)

# identify mean and get closest N animals
mo6$dist <- abs(mo6$Protexpr - cal[cal$Age == 6,]$Protexpr)
mo18$dist <- abs(mo18$Protexpr - cal[cal$Age == 18,]$Protexpr)
mo6 <- arrange(mo6, dist)
mo18 <- arrange(mo18, dist)

# combine and sort
select <- rbind(mo6[1:Num,], mo18[1:Num,])
select <- arrange(select, Age, Mouse.ID)
print(select)

# is selected?
samples$is_select <- NA
samples$is_select <- samples$Mouse.ID %in% select$Mouse.ID

# Plot all 6 and 18 months
ggplot(samples, aes(x = RNAexpr, y = Protexpr, colour = as.factor(as.character(Age)))) +
  geom_point() +
  geom_smooth( method = lm, se = FALSE) +
  labs( title = "Females") +
  scale_color_aaas()
