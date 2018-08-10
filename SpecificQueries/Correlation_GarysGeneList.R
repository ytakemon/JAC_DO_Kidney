# R/3.4.4

# load libraries
library(tidyverse)
library(corrplot)

# set up directories
basedir <- "/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/"
setwd(basedir)

# Load data
load("./RNAseq_data/DO188b_kidney.RData")
rm(G, Glist, N, covar, genoprobs, raw.mrna, raw.protein, snps)
GeneList <- read.csv("./GaryList.csv", stringsAsFactors = FALSE)
names(GeneList) <- c("Symbol", "Category")

# Check list
any((GeneList$Symbol %in% annot.mrna$symbol) == FALSE)
filter(GeneList, (Symbol %in% annot.mrna$symbol) == FALSE)
#   Symbol   Category
# 1    Lrp ProxTubule ---> Lrp2
# 2  Lamc4       Glom ---> ??? No lamc4 gene... remove
# 3 Ndfa11     Oxphos ---> Ndufa11

# Fix list
GeneList[GeneList$Symbol == "Lrp",]$Symbol <- "Lrp2"
GeneList[GeneList$Symbol == "Ndfa11",]$Symbol <- "Ndufa11"
GeneList <- GeneList[GeneList$Symbol != "Lamc4",]

# Check list again
any((GeneList$Symbol %in% annot.mrna$symbol) == FALSE)

# Annotate with ENSEMBL ID
GeneList$gene_id <- NA
GeneList$prot_id <- NA
for(i in 1:nrow(GeneList)){
  GeneList$gene_id[i] <- annot.mrna[annot.mrna$symbol == GeneList$Symbol[i],]$id

  if(any(annot.protein$symbol == GeneList$Symbol[i]) == FALSE){
    GeneList$prot_id[i] <- NA
  } else{
    GeneList$prot_id[i] <- annot.protein[annot.protein$symbol == GeneList$Symbol[i],]$id
  }
}

# get sublist
for(type in c("mrna", "protein")){

  # Determine datasets
  if(type == "mrna"){
    List <- GeneList
    Expr <- expr.mrna[,List$gene_id]
  } else if(type == "protein"){
    List <- GeneList[!is.na(GeneList$prot_id),]
    Expr <- expr.protein[,List$prot_id]
  }

  # Get ages
  for(age in c("all","6mo","12mo","18mo")){

    # Determine age and get pearson correlation
    if(age == "all"){
      ExprCor <- cor(Expr)
      rownames(ExprCor) <- colnames(ExprCor) <- List$Symbol
    } else if(age == "6mo"){
      ExprCor <- cor(Expr[annot.samples[annot.samples$Age == 6,]$Mouse.ID,])
      rownames(ExprCor) <- colnames(ExprCor) <- List$Symbol
    } else if(age == "12mo"){
      ExprCor <- cor(Expr[annot.samples[annot.samples$Age == 12,]$Mouse.ID,])
      rownames(ExprCor) <- colnames(ExprCor) <- List$Symbol
    } else if(age == "18mo"){
      ExprCor <- cor(Expr[annot.samples[annot.samples$Age == 18,]$Mouse.ID,])
      rownames(ExprCor) <- colnames(ExprCor) <- List$Symbol
    }

    filename <- paste0("./Plot/GaryGeneList_",type,"_",age,"_age.pdf")
    title <- paste0(type," (RankZ) pearson correlation at age: ", age)
    # Plot
    set.seed(123)
    pdf(filename, height = 8, width = 8)
    corrplot(ExprCor,
        method = "circle",
        order = "hclust",
        title = title,
        mar = c(0,0,2,0),
        tl.col = "black")
    dev.off()
  }
}





# get sublist
for(type in c("mrna", "protein")){

  # Determine datasets
  if(type == "mrna"){
    List <- GeneList
    Expr <- expr.mrna[,List$gene_id]
  } else if(type == "protein"){
    List <- GeneList[!is.na(GeneList$prot_id),]
    Expr <- expr.protein[,List$prot_id]
  }

  # Get ages
  for(age in c("all","6mo","12mo","18mo")){

    # Determine age and get pearson correlation
    if(age == "all"){
      # Get correlation
      set.seed(123)
      ExprCor <- cor(Expr)
      rownames(ExprCor) <- colnames(ExprCor) <- List$Symbol

      # get all age order
      symbol <- rownames(ExprCor)
      order <- symbol[corrMatOrder(ExprCor, order = "hclust")]

      # reorder
      ExprCor <- ExprCor[order,order]

    } else if(age == "6mo"){
      # Get correlation
      ExprCor <- cor(Expr[annot.samples[annot.samples$Age == 6,]$Mouse.ID,])
      rownames(ExprCor) <- colnames(ExprCor) <- List$Symbol

      # reorder
      ExprCor <- ExprCor[order,order]

    } else if(age == "12mo"){
      # Get correlation
      ExprCor <- cor(Expr[annot.samples[annot.samples$Age == 12,]$Mouse.ID,])
      rownames(ExprCor) <- colnames(ExprCor) <- List$Symbol

      # reorder
      ExprCor <- ExprCor[order,order]

    } else if(age == "18mo"){
      # Get correlation
      ExprCor <- cor(Expr[annot.samples[annot.samples$Age == 18,]$Mouse.ID,])
      rownames(ExprCor) <- colnames(ExprCor) <- List$Symbol

      # reorder
      ExprCor <- ExprCor[order,order]

    }

    filename <- paste0("./Plot/GaryGeneList_",type,"_",age,"_age.pdf")
    title <- paste0(type," (RankZ) pearson correlation at age: ", age)

    # Plot
    set.seed(123)
    pdf(filename, height = 8, width = 8)
    corrplot(ExprCor,
        method = "circle",
        title = title,
        mar = c(0,0,2,0),
        tl.col = "black")
    dev.off()
  }
}
