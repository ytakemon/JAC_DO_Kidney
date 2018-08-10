library(tidyverse)
library(readxl)
library(corrplot)
library(ggsci)
# change global options
options(tibble.width = Inf)

# Get data
load("~/Dropbox/TheAgingKidneyData/RData/DO188b_kidney_noprobs.RData")
GaryList <- read_excel("~/Dropbox/JAC_DO_Kidney_Aging/Other analyses not in manuscript/Pathway/Correlations/GeneList_fromGary.xlsx")

# remove duplicates
GaryList <- GaryList[!duplicated(GaryList$GarysGeneList),]
dim(GaryList) # total is 75

# Subset to RNA
GaryList_RNA <- GaryList_RNA <- annot.mrna[annot.mrna$symbol %in% GaryList$GarysGeneList,]
dim(GaryList_RNA) # 71

# Adjust Gary's GaryList
GaryList_adj <- GaryList[GaryList$GarysGeneList %in% GaryList_RNA$symbol,]
dim(GaryList_adj) # 71

# Check out protein protein
GaryList_protein <- GaryList_adj[GaryList_adj$GarysGeneList %in% annot.protein$symbol,]
dim(GaryList_protein) # 52

# annotate protein exists
for(i in 1:nrow(GaryList_adj)){
  if(any(GaryList_protein$GarysGeneList == GaryList_adj$GarysGeneList[i])){
    GaryList_adj$ProteinData[i] <- TRUE
  } else {
    GaryList_adj$ProteinData[i] <- FALSE
  }
}

# grab expression annotation and add category
# m
sub_annot_m <- annot.mrna[annot.mrna$symbol %in% GaryList_adj$GarysGeneList,]
sub_annot_m$Category <- NA
for(i in 1:nrow(sub_annot_m)){
  sub_annot_m$Category[i] <- GaryList_adj[GaryList_adj$GarysGeneList == sub_annot_m$symbol[i],]$Category
}
sub_annot_m <- arrange(sub_annot_m, Category)
# p
sub_annot_p <- annot.protein[annot.protein$symbol %in% GaryList_adj[GaryList_adj$ProteinData == TRUE,]$GarysGeneList,]
sub_annot_p$Category <- NA
for(i in 1:nrow(sub_annot_p)){
  sub_annot_p$Category[i] <- GaryList_adj[GaryList_adj$GarysGeneList == sub_annot_p$symbol[i],]$Category
}
sub_annot_p <- arrange(sub_annot_p, Category)

# subset and reannotate with symbol
# m
expr_m <- expr.mrna[, sub_annot_m$id]
colnames(expr_m) <- sub_annot_m$symbol
# p
expr_p <- expr.protein[, sub_annot_p$id]
colnames(expr_p) <- sub_annot_p$symbol

for(i in c("all", 6, 12, 18)){
  if(i == "all"){
    # All ages:
    # m:
    pdf("~/Desktop/GarySelect_mrna_ALLage.pdf", width = 12, height = 12)
    corrplot(cor(expr_m), method = "circle", tl.col = "black")
    dev.off()
    # p:
    pdf("~/Desktop/GarySelect_protein_ALLage.pdf", width = 12, height = 12)
    corrplot(cor(expr_p, use = "complete.obs" ), method = "circle", tl.col = "black")
    dev.off()
  } else if(i == 6){
    # subset
    sub_expr_m <- expr_m[annot.samples$Age == 6,]
    sub_expr_p <- expr_p[annot.samples$Age == 6,]

    # m:
    pdf("~/Desktop/GarySelect_mrna_6MOage.pdf", width = 12, height = 12)
    corrplot(cor(sub_expr_m), method = "circle", tl.col = "black")
    dev.off()
    # p:
    pdf("~/Desktop/GarySelect_protein_6MOage.pdf", width = 12, height = 12)
    corrplot(cor(sub_expr_p, use = "complete.obs" ), method = "circle", tl.col = "black")
    dev.off()
  } else if(i == 12){
    # subset
    sub_expr_m <- expr_m[annot.samples$Age == 12,]
    sub_expr_p <- expr_p[annot.samples$Age == 12,]

    # m:
    pdf("~/Desktop/GarySelect_mrna_12MOage.pdf", width = 12, height = 12)
    corrplot(cor(sub_expr_m), method = "circle", tl.col = "black")
    dev.off()
    # p:
    pdf("~/Desktop/GarySelect_protein_12MOage.pdf", width = 12, height = 12)
    corrplot(cor(sub_expr_p, use = "complete.obs" ), method = "circle", tl.col = "black")
    dev.off()
  } else if(i == 18){
    # subset
    sub_expr_m <- expr_m[annot.samples$Age == 18,]
    sub_expr_p <- expr_p[annot.samples$Age == 18,]

    # m:
    pdf("~/Desktop/GarySelect_mrna_18MOage.pdf", width = 12, height = 12)
    corrplot(cor(sub_expr_m), method = "circle", tl.col = "black")
    dev.off()
    # p:
    pdf("~/Desktop/GarySelect_protein_18MOage.pdf", width = 12, height = 12)
    corrplot(cor(sub_expr_p, use = "complete.obs" ), method = "circle", tl.col = "black")
    dev.off()
  }
}

# Pull out genes for mRNA comparison

get.gene <- function(gene.symbol){
  # pull covariates and convert to factors
  # note Age is numeric and fAge is factor
  Mouse.ID <- rownames(annot.samples)
  Sex <- factor(as.character(model.matrix(~annot.samples[,"Sex"])[,2]+1), levels=c("1","2"), labels=c("F","M"))
  Age <- annot.samples$Age
  fAge <- factor(as.character(Age),levels=c("6","12","18"))

  #pull gene identifiers from mrna and protein annotations
  mrna.id <- annot.mrna[annot.mrna$symbol==gene.symbol,"id"]
  protein.id <- annot.protein[annot.protein$symbol==gene.symbol,"id"]

  # strat building the tiblle to hold data
  my.data <- data_frame(Mouse.ID=Mouse.ID,
                        Sex=Sex,
                        Age=Age,
                        fAge=fAge)
  # add RNA and Protein columns to the tibble - if it exists
  if(length(mrna.id)>0){
    my.data <- mutate(my.data, RNA=expr.mrna[,mrna.id])
  }
  if(length(protein.id)>0){
    my.data <- mutate(my.data, Protein=expr.protein[,protein.id])
  }
  # return value is a tibble with covariates, RNA and Protein data
  my.data
}

combine_genes <- function(gene1, gene2, type){
  # grab both df
  df1 <- get.gene(gene1)
  df2 <- get.gene(gene2)
  # combine and rename columns
  if(type == "mRNA"){
    df1 <- df1[,c("Mouse.ID","Sex","Age","fAge","RNA")]
    df2 <- df2[,c("Mouse.ID","Sex","Age","fAge","RNA")]
    df <- cbind(df1, df2$RNA)
    colnames(df) <- c("Mouse.ID","Sex","Age","fAge","Gene1","Gene2")
  } else if(type == "Protein"){
    if(length(colnames(df1)) == 6 & length(colnames(df2)) == 6){
      df1 <- df1[,c("Mouse.ID","Sex","Age","fAge","Protein")]
      df2 <- df2[,c("Mouse.ID","Sex","Age","fAge","Protein")]
      df <- cbind(df1, df2$Protein)
      colnames(df) <- c("Mouse.ID","Sex","Age","fAge","Gene1","Gene2")
    } else if(length(colnames(df1)) < 6 & length(colnames(df2)) == 6){
      df1 <- df1[,c("Mouse.ID","Sex","Age","fAge")]
      df2 <- df2[,c("Mouse.ID","Sex","Age","fAge","Protein")]
      df <- cbind(df1, df2$Protein)
      colnames(df) <- c("Mouse.ID","Sex","Age","fAge", paste0(gene2,"_Protein"))
    } else if(length(colnames(df1)) == 6 & length(colnames(df2)) < 6){
      df1 <- df1[,c("Mouse.ID","Sex","Age","fAge","Protein")]
      colnames(df) <- c("Mouse.ID","Sex","Age","fAge", paste0(gene1,"_Protein"))
    }
  }
  # return values of RNA and Protei from both genes
  df
}

createPlot <- function(gene1, gene2, type){
  df <- combine_genes(gene1, gene2, type)

  # check point
  if(ncol(df) < 6){
    # spit out df and stop
    df
    stop("Need two genes for comparison")
  }

  df <- df %>%
    mutate(Fitted = fitted(lm(Gene2 ~ Gene1 + Sex + fAge + Sex:fAge, data=.)))
  df_summary <- df %>%
    group_by(Sex, fAge) %>%
    summarise(mGene1 = mean(Gene1), sdGene1 = sd(Gene1),
             mGene2 = mean(Gene2), sdGene2 = sd(Gene2),
             N = n()) %>%
    mutate(seGene1 = sdGene1/sqrt(N), seGene2 = sdGene2/sqrt(N))
  df_plot <- ggplot(df_summary, aes(x = mGene1, y = mGene2, colour = fAge, shape = Sex), size = 4) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = mGene2 - seGene2, ymax = mGene2 + seGene2), size = 1)+
    geom_errorbarh(aes(xmin = mGene1 - seGene1, xmax = mGene1 + seGene1), size = 1)+
    geom_point(data = df, aes(x = Gene1, y = Gene2))+
    geom_line(data = df, aes(x = Gene1, y = Fitted))+
    xlab(paste(gene1, type))+
    ylab(paste(gene2, type))+
    theme_bw() +
    scale_colour_aaas()
  # return plot
  df_plot
}

# same direction
# Proximal tub - Proximal tub
# Ass1 and Acsm2
plot1 <- createPlot("Ass1","Acsm2", type = "mRNA")
pdf("~/Desktop/temp/Ass1_Acsm2_mRNA_comparison.pdf", height = 5, width = 6)
plot1
dev.off()
# Prox - Oxphos
plot2 <- createPlot("Ass1","Pink1", type = "mRNA")
pdf("~/Desktop/temp/Ass1_Pink1_mRNA_comparison.pdf", height = 5, width = 6)
plot2
dev.off()
# Prox - Glom
plot3 <- createPlot("Ass1","Podxl", type = "mRNA")
pdf("~/Desktop/temp/Ass1_Podxl_mRNA_comparison.pdf", height = 5, width = 6)
plot3
dev.off()
# Oxphos - Glom
plot4 <- createPlot("Pink1","Podxl", type = "mRNA")
pdf("~/Desktop/temp/Pink1_Podxl_mRNA_comparison.pdf", height = 5, width = 6)
plot4
dev.off()
