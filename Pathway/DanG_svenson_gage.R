################################################################################
# Search for significantly altered functional categories of genes in the liver 
# GE data.
# Daniel Gatti
# dan.gatti@jax.org
# Oct. 22, 2016
################################################################################
options(stringsAsFactors = FALSE)
library(gage)
library(AnnotationHub)
library(biomaRt)
library(KEGGREST)
library(org.Mm.eg.db)
library(doParallel)

base_dir = "/projects/churchill-lab/projects/Svenson_DO850/Gatti_liver_physiol/"
setwd(base_dir)

out_dir = paste0(base_dir, "results/gage/KEGG")

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

# Register parallel back end.
registerDoParallel(cores = 19)

# Set up covariates.
covar = model.matrix(~sex + diet + gen + litter, data = dataset.expr$samples)[,-1]
stopifnot(colnames(expr.mrna) == rownames(covar))

# Sex & Diet.
# Regress out generation & litter.
expr.res = expr.mrna
for(i in 1:nrow(expr.mrna)) {

  mod = lm(expr.mrna[i,] ~ covar[,3:8])
  expr.res[i,] = residuals(mod)

} # for(i)

obj = gage(expr.res, gsets = kegg.paths, ref = which(covar[,"diethf"] == 0), 
           samp = which(covar[,"diethf"] == 1), compare = "unpaired",
           set.size = c(5, 200), use.fold = FALSE)


