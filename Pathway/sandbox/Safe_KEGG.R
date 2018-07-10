# SAFE - KEGG analysis adapted from Dan Gatti's svenson_SAFE_KEGG.R
library(dplyr)
library(safe)
library(GOstats)
library(Rgraphviz)
library(AnnotationHub)
library(biomaRt)
library(KEGGREST)
library(org.Mm.eg.db)
library(doParallel)

# set up directories
basedir <- "/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/"
setwd(basedir)

outdir <- paste0(basedir,"Pathways/")

# Get Ensembl genes
hub <- AnnotationHub()
hub <- query(hub, c("ensembl", "mus musculus"))
ensembl <- hub[[names(hub)[hub$title == "Mus_musculus.GRCm38.82.gtf"]]]

# load in rankZ normalied expression data
load("./RNAseq_data/DO188b_kidney.RData")
rm(genoprobs, G, Glist, N, raw.protein, raw.mrna, snps)

# Subset Ensembl data
ensembl <- ensembl[ensembl$type == "gene"]
ensembl <- ensembl[ensembl$gene_id %in% colnames(expr.mrna)]
mrna <- expr.mrna[,ensembl$gene_id]
stopifnot(colnames(mrna) == ensembl$gene_id)

# Get Entrez IDs for the Ensembl IDs in our data. KEGG uses Entrez Gene IDs.
ens2entrez <- AnnotationDbi::select(org.Mm.eg.db, key = ensembl$gene_id,
             column = c("ENSEMBL", "ENTREZID","SYMBOL"), keytype = "ENSEMBL")
colnames(ens2entrez) <- tolower(colnames(ens2entrez))

# Get KEGG terms for the ensembl genes.
kegg <- keggLink("pathway", "mmu")
kegg <- data.frame(path = kegg, entrezid = names(kegg))

# Get the KEGG pathway names for mouse.
kegg.desc <- keggList("pathway", "mmu")
kegg.desc <- data.frame(desc = kegg.desc, path = names(kegg.desc))

# merge
kegg <- merge(kegg, kegg.desc, by = "path")
kegg <- kegg %>% mutate(path      = sub("^path:", "", kegg$path),
                        entrezid  = sub("^mmu:", "", kegg$entrezid),
                        desc      = sub(" - Mus musculus \\(mouse\\)$", "", kegg$desc))
kegg <- merge(kegg, ens2entrez, by = "entrezid")

# Convert this ensembl ID based list to a KEGG pathway based list.
kegg.paths <- split(kegg$ensembl, kegg$path)

# Calculate the sex group mean for each KEGG pathway
sex <- factor(annot.samples$Sex)
names(sex) <- rownames(annot.samples)
mrna <- t(expr.mrna)
stopifnot(names(sex) == colnames(mrna))

kegg_means <- matrix(data = 0,
                     nrow = length(kegg.paths),
                     ncol = 2)
rownames(kegg_means) <- names(kegg.paths)
colnames(kegg_means) <- levels(sex)

for (i in seq_along(kegg.paths)){
  g <- t(mrna[kegg.paths[[i]],])
  g <- data.frame(g)

  if(nrow(g) == 1){
    g <- t(g)
    colnames(g) <- kegg.paths[[i]]
    rownames(g) <- colnames(mrna)
    g <- as.data.frame(g)
    g <- split(g, t(sex))
    Fem <- mean(g$F[,1])
    Mal <- mean(g$M[,1])
    kegg_means[i,"F"] <- Fem
    kegg_means[i,"M"] <- Mal
  } else {
    g <- split(g, t(sex))
    kegg_means[i,] = colMeans(sapply(g, colMeans, na.rm = TRUE))
  }
}

# Register parallel back end.

# Custom function fors sex * age interaction.
local.inter.model = function(X.mat, y.vec, ...) {

  sex  = as.numeric(substring(y.vec, 1, 1))
  diet = as.numeric(substring(y.vec, 3, 3))
  local.covar = model.matrix(~ Sex + Age + Sex*Age, data = annot.samples)[,-1]

  return(
    function(data, vector = local.covar, ...) {

      # Sex + Age model
      qry = qr(local.covar[,1:3])
      resid.null = apply(data, 1, function(z) { qr.resid(qry, z) })
      # Sex + Age + Sex*Age model
      qry = qr(local.covar)
      resid.full = apply(data, 1, function(z) { qr.resid(qry, z) })
      return(-ncol(data) * log(colSums(resid.full^2) / colSums(resid.null^2)))
    }
  )
} # local.inter.model()

safe.fxn = function(expr, kegg.list, y.vec, file.prefix, model = "default",
           Z.mat = NULL) {

  C.mat = getCmatrix(kegg.list, present.genes = rownames(expr), min.size = 5,
          max.size = 600)

  # Subset the expression again to match the C.matrix.
  expr = expr[C.mat$row.names,]
  print(paste(nrow(expr), "genes &", ncol(C.mat$C.mat.csr), "categories"))
  stopifnot(rownames(expr) == C.mat$row.names)

  results = safe(X.mat = expr, y.vec = y.vec, C.mat = C.mat, local = model,
            alpha = 1.0, error = "FDR.BH", method = "bootstrap", Pi.mat = 10, parallel = TRUE)

  saveRDS(results, file = paste0(file.prefix, ".rds"))

  pv = results@global.pval
  nm = kegg$desc[match(names(pv), kegg$path)]

  pv = pv[!is.na(nm)]
  nm = nm[!is.na(nm)]

  table = data.frame(KEGG = names(pv), Name = nm, pv)
  saveRDS(table, file = paste0(file.prefix, "_table.rds"))

} # safe.fxn()


# Set up covariates
covariates <- model.matrix(~ Sex + Age + Generation, data = annot.samples)[,-1]
stopifnot(colnames(mrna) == rownames(covariates))

# Sex only one-way ANOVA,
# regress out age and gen first
expr.res <- mrna
for(i in 1:nrow(mrna)) {
  print(i)
  mod = lm(mrna[i,] ~ covariates[,2:ncol(covariates)])
  expr.res[i,] = residuals(mod)
}

registerDoParallel(cores = 10)
y.vec <- covariates[,"SexM"]
safe.fxn( expr = expr.res,
          kegg.list = kegg.paths,
          y.vec = y.vec,
          file.prefix = paste0(outdir, "SexOnly"),
          model = "default")

# Age only one-way ANOVA,
# regress out sex and gen first
expr.res <- mrna
for(i in 1:nrow(mrna)) {
  print(i)
  mod = lm(mrna[i,] ~ covariates[,c(1,3:ncol(covariates))])
  expr.res[i,] = residuals(mod)
}

registerDoParallel(cores = 10)
y.vec <- covariates[,"Age"]
safe.fxn( expr = expr.res,
          kegg.list = kegg.paths,
          y.vec = y.vec,
          file.prefix = paste0(outdir, "AgeOnly"),
          model = "default")

# Sex * Age interaction ANOVA,
# regress out generation.
expr.res <- mrna
for(i in 1:nrow(mrna)) {
  print(i)
  mod = lm(mrna[i,] ~ covariates[,3:ncol(covariates)])
  expr.res[i,] = residuals(mod)
}

registerDoParallel(cores = 10)
y.vec <- apply(covariates[,c("Age","SexM")], 1, paste0, collapse = ":")
safe.fxn( expr = expr.res,
          kegg.list = kegg.paths,
          y.vec = y.vec,
          file.prefix = paste0(outdir, "SexAgeInteraction"),
          model = "inter.model")

# Collect all tables generated
files <- list.files(path = "./Pathways", pattern = "_table.rds$", full.names = TRUE)
for (file in files){
  table <- readRDS(file)
  obj <- readRDS(sub("_table", "", file))
  err <- as.matrix(obj@global.error)
  table <- merge(table, err, by = "row.names")
  rownames(table) <- table[,1]
  colnames(table)[5] <- "FDR_BH"
  table <- table[,-1]
  table[,4] <- format(table[,4], digit = 4)
  table <- table[order(table[,3]),]
  write.table(table, file = sub("rds$","txt", file), sep = "\t", quote = FALSE, row.names = FALSE)
}

















##########
# Harvest the tables.
setwd(out_dir)

files = dir(pattern = "_table.rds$")
for(f in files){
  tbl = readRDS(f)
  obj = readRDS(sub("_table", "", f))
  err = as.matrix(obj@global.error)
  tbl = merge(tbl, err, by = "row.names")
  rownames(tbl) = tbl[,1]
  colnames(tbl)[5] = "FDR.BH"
  tbl = tbl[,-1]
  tbl[,4] = format(tbl[,4], digit = 4)
  tbl = tbl[order(tbl[,3]),]
  write.table(tbl, file = sub("rds$", "txt", f), sep = "\t", quote = F,
              row.names = F)
} # for(f)


##########
# Add the mean expression in each sex/diet group.
gm = readRDS("safe_KEGG_cat_gene_means.rds")

files = dir(pattern = "_GOBP_.+_table.txt$")

for(f in files){

  tbl = read.delim(f)
  rownames(tbl) = tbl[,1]
  tbl = merge(tbl, gm, by = "row.names")
  tbl = tbl[,-1]
  tbl = tbl[order(tbl$pv),]
  write.table(tbl, file = f, sep = "\t", quote = F, row.names = F)

} # for(f)




##########
# Make SAFE plots for categories with p < 0.05.
obj.files = dir(pattern = ".rds$")
obj.files = obj.files[-grep("table|means", obj.files)]

for(i in 1:length(obj.files)) {

  tbl.file = sub("\\.rds$", "_table.rds", obj.files[i])
  obj = readRDS(obj.files[i])
  tbl = readRDS(tbl.file)
  tbl = tbl[tbl[,3] <= 0.05,]
  tbl = tbl[order(tbl[,3]),]

  pdf(sub("rds", "pdf", obj.files[i]), width = 8, height = 6)
  for(j in 1:nrow(tbl)) {
    safeplot(obj, rownames(tbl)[j])
  } # for(j)
  dev.off()

} # for(i)
