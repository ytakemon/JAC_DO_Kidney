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
Num <- 3

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

# identify mean
min6 <- which.min(abs(mo6$expr - cal[cal$Age == 6,]$expr))
min18 <- which.min(abs(mo18$expr - cal[cal$Age == 18,]$expr))

if( Num == 1){
  Num <- 0
  min6 <- min6 + Num
  min18 <- min18 + Num
} else if( Num == 2){
  Num <- c(0,1)
  min6 <- min6 + Num
  min18 <- min18 + Num
} else if( Num == 3){
  Num <- c(-1, 0, 1)
  min6 <- min6 + Num
  min18 <- min18 + Num
} else if( Num == 4){
  Num <- c(-1, 0, 1, 2)
  min6 <- min6 + Num
  min18 <- min18 + Num
} else if( Num == 5){
  Num <- c(-2, -1, 0, 1, 2)
  min6 <- min6 + Num
  min18 <- min18 + Num
} else if( Num == 6){
  Num <- c(-2, -1, 0, 1, 2, 3)
  min6 <- min6 + Num
  min18 <- min18 + Num
} else if( Num == 7){
  Num <- c(-3, -2, -1, 0, 1, 2, 3)
  min6 <- min6 + Num
  min18 <- min18 + Num
}
min6 <- mo6[min6,]
min18 <- mo18[min18,]
select <- rbind(min6, min18)

# Print results
print(select)
