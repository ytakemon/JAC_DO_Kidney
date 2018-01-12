# Identify animals that have the average given expression values for the given group.
library(dplyr)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")
load("./shiny_annotation.RData")
erk <- read.delim("./Phenotype/phenotypes/JAC_WB_kidney_ERK.txt")
akt <- read.delim("./Phenotype/phenotypes/JAC_WB_kidney_AKT.txt")

# Selection quality:
Gene_name <- "Lrp2"
Data <- "protein"
Sex <- "M"
Age <- c(6, 18)
Num <- 5
Type <- "omics" #/ "WB_Erk1" / "WB_Akt1"
Query <- "max" # "min"
erk1 <- erk
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

if ( Type == "omics"){
  gene <- other.ids(Gene_name, Data)
  # Get expression data
  if(Data == "mRNA"){
    samples$expr <- expr.mrna[, gene$id]
  } else if (Data == "protein"){
    samples$expr <- expr.protein[, gene$protein_id]
  }

  # Subset selction
  samples <- samples[ samples$Age %in% Age,]
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

} else if (Type == "WB_Erk1" & Query == "max"){
  erk1 <- erk1[complete.cases(erk1$Phospho_ERK1),]
  erk1 <- erk1[erk1$ID %in% samples$Mouse.ID,]
  samples <- samples[samples$Mouse.ID %in% erk1$ID,]
  samples$expr <- erk1$Phospho_ERK1
  samples <- arrange(samples, desc(expr))
  samples <- samples[(samples$Age == 18) & (samples$Sex == "M"),]
  head(samples)
}
