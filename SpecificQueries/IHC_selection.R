# Identify animals that have the average given expression values for the given group.
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")
load("./shiny_annotation.RData")

# Selection quality:
Gene_name <- "Lrp2"
Data <- "protein"
Sex <- "M"
Age <- c(6, 18)
Num <- 5

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
if(Data == "mRNA"){
  annot.samples$expr <- expr.mrna[, gene$id]
} else if (Data == "protein"){
  annot.samples$expr <- expr.protein[, gene$protein_id]
}

# Subset selction
samples <- annot.samples[ annot.samples$Age %in% Age,]
samples <- samples[ samples$Sex %in% Sex,]
cal <- aggregate(expr ~ Age, data = samples, mean)
mo6 <- samples[ samples$Age == 6,] %>% arrange(expr)
mo18 <- samples[ samples$Age == 18,] %>% arrange(expr)

# identify mean and get closest N animals
mo6$dist <- abs(mo6$expr - cal[cal$Age == 6,]$expr)
mo18$dist <- abs(mo18$expr - cal[cal$Age == 18,]$expr)
mo6 <- arrange(mo6, dist)
mo18 <- arrange(mo18, dist)

# combine and sort
select <- rbind(mo6[1:Num,], mo18[1:Num,])
select <- arrange(select, Age, Mouse.ID)
print(select)
