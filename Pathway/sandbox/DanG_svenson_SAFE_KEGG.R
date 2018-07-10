################################################################################
# Search for significantly altered functional categories of genes in the liver 
# GE data.
# Daniel Gatti
# dan.gatti@jax.org
# Oct. 22, 2016
################################################################################
library(safe)
library(GOstats)
library(Rgraphviz)
library(AnnotationHub)
library(biomaRt)
library(KEGGREST)
library(org.Mm.eg.db)
library(doParallel)

base_dir = "/projects/churchill-lab/projects/Svenson_DO850/Gatti_liver_physiol/"
setwd(base_dir)

out_dir = paste0(base_dir, "results/SAFE/KEGG/")

# Get Ensembl genes.
hub = AnnotationHub()
hub = query(hub, c("ensembl", "mus musculus"))
ensembl = hub[[names(hub)[hub$title == "Mus_musculus.GRCm38.82.gtf"]]]

# Load in the rankZ normalized expression data.
load("data/Svenson_DO850_for_eQTL_veiwer_v0.Rdata")
rm(genoprobs, K, dataset.expr.do384, dataset.expr.petr)

# Transpose the expression data for SAFE.
# Using the normalized expression data, not the rankZ data.
expr.mrna = t(dataset.expr$norm)

# Subset the Ensembl data to only include the genes in the DO478 data set.
ensembl = ensembl[ensembl$type == "gene"]
ensembl = ensembl[ensembl$gene_id %in% rownames(expr.mrna)]
expr.mrna = expr.mrna[ensembl$gene_id,]
stopifnot(rownames(expr.mrna) == ensembl$gene_id)

# Get Entrez IDs for the Ensembl IDs in our data. KEGG uses Entrez Gene IDs.
ens2entrez = AnnotationDbi::select(org.Mm.eg.db, key = ensembl$gene_id, 
             column = c("ENSEMBL", "ENTREZID","SYMBOL"), keytype = "ENSEMBL")
colnames(ens2entrez) = tolower(colnames(ens2entrez))

# Get KEGG terms for the ensembl genes.
kegg = keggLink("pathway", "mmu")
kegg = data.frame(path = kegg, entrezid = names(kegg))

# Get the KEGG pathway names for mouse.
kegg.desc = keggList("pathway", "mmu")
kegg.desc = data.frame(desc = kegg.desc, path = names(kegg.desc))

kegg = merge(kegg, kegg.desc, by = "path")
kegg$path     = sub("^path:", "", kegg$path)
kegg$entrezid = sub("^mmu:", "", kegg$entrezid)
kegg$desc     = sub(" - Mus musculus \\(mouse\\)$", "", kegg$desc)

kegg = merge(kegg, ens2entrez, by = "entrezid")

# Convert this ensembl ID based list to a KEGG pathway based list.
kegg.paths = split(kegg$ensembl, kegg$path)

length(kegg.paths)

# Calculate the sex/diet group means for each KEGG pathway.
sex_diet = factor(paste0(dataset.expr$samples$sex, dataset.expr$samples$diet))
names(sex_diet) = dataset.expr$samples$sample
stopifnot(names(sex_diet) == colnames(expr.mrna))

kegg_means = matrix(0, nrow = length(kegg.paths), ncol = 4, dimnames = 
             list(names(kegg.paths), levels(sex_diet)))

for(i in seq_along(kegg.paths)) {

  g = expr.mrna[kegg.paths[[i]],]
  g = data.frame(t(g))
  g = split(g, sex_diet)
  kegg_means[i,] = colMeans(sapply(g, colMeans, na.rm = T))

} # for(i)

saveRDS(kegg_means, file = paste0(out_dir, "safe_KEGG_cat_gene_means.rds"))

# Custom function for sex X diet interaction.
local.inter.model = function(X.mat, y.vec, ...) {

  sex  = as.numeric(substring(y.vec, 1, 1))
  diet = as.numeric(substring(y.vec, 3, 3))
  local.covar = model.matrix(~ sex + diet + sex*diet)

  return(
    function(data, vector = local.covar, ...) {

      # Sex + Diet model
      qry = qr(local.covar[,1:3])
      resid.null = apply(data, 1, function(z) { qr.resid(qry, z) })
      # Sex + Diet + Sex*Diet model
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
            alpha = 1.0, error = "FDR.BH", method = "bootstrap", Pi.mat = 10000, parallel = TRUE)

  saveRDS(results, file = paste0(file.prefix, ".rds"))

  pv = results@global.pval
  nm = kegg$desc[match(names(pv), kegg$path)]

  pv = pv[!is.na(nm)]
  nm = nm[!is.na(nm)]

  table = data.frame(KEGG = names(pv), Name = nm, pv)
  saveRDS(table, file = paste0(file.prefix, "_table.rds"))

} # safe.fxn()


# Register parallele back end.
registerDoParallel(cores = 19)

# Set up covariates.
covar = model.matrix(~sex + diet + gen + litter, data = dataset.expr$samples)[,-1]
stopifnot(colnames(expr.mrna) == rownames(covar))

# Sex & Diet one-way ANOVA.
expr.res = expr.mrna
y.vec = apply(covar[,c("sexM", "diethf")], 1, paste, collapse = ":")
safe.fxn(expr = expr.res, kegg.list = kegg.paths, y.vec = y.vec,
         file.prefix = paste0(out_dir, "DO478_safe_KEGG_sex_diet"))

# Sex only, regress out diet and gen first.
expr.res = expr.mrna
for(i in 1:nrow(expr.mrna)) {

  mod = lm(expr.mrna[i,] ~ covar[,2:8])
  expr.res[i,] = residuals(mod)

} # for(i)

y.vec = covar[,"sexM"]
safe.fxn(expr = expr.res, kegg.list = kegg.paths, y.vec = y.vec,
         file.prefix = paste0(out_dir, "DO478_safe_KEGG_sex"))

# Diet only, regress out sex and gen first.
expr.res = expr.mrna
for(i in 1:nrow(expr.mrna)) {

  mod = lm(expr.mrna[i,] ~ covar[,c(1,3:8)])
  expr.res[i,] = residuals(mod)

} # for(i)

y.vec = covar[,"diethf"]
safe.fxn(expr = expr.res, kegg.list = kegg.paths, y.vec = y.vec,
         file.prefix = paste0(out_dir, "DO478_safe_KEGG_diet"))

# Sex X Diet ANOVA.
# Regress out generation.
expr.res = expr.mrna
for(i in 1:nrow(expr.mrna)) {

  mod = lm(expr.mrna[i,] ~ covar[,3:8])
  expr.res[i,] = residuals(mod)

} # for(i)
y.vec = apply(covar[,c("sexM", "diethf")], 1, paste, collapse = ":")

safe.fxn(expr = expr.res, kegg.list = kegg.paths, y.vec = y.vec,
         file.prefix = paste0(out_dir, "DO478_safe_KEGG_sexXdiet"),
         model = "inter.model")


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




