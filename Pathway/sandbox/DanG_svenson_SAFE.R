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
library(GO.db)
library(doParallel)

base_dir = "/projects/churchill-lab/projects/Svenson_DO850/Gatti_liver_physiol/"
setwd(base_dir)

out_dir = paste0(base_dir, "results/SAFE/")

# Get Ensembl genes.
hub = AnnotationHub()
hub = query(hub, c("ensembl", "mus musculus"))
ensembl = hub[[names(hub)[hub$title == "Mus_musculus.GRCm38.82.gtf"]]]

# Load in the rankZ normalized expression data.
load("data/Svenson_DO850_for_eQTL_veiwer_v0.Rdata")
rm(genoprobs, K, dataset.expr.do384, dataset.expr.petr)

# Transpose the expression data for SAFE.
expr.mrna = t(dataset.expr$expr)

# Subset the Ensembl data to only include the genes in the DO478 data set.
ensembl = ensembl[ensembl$type == "gene"]
ensembl = ensembl[ensembl$gene_id %in% rownames(expr.mrna)]
expr.mrna = expr.mrna[ensembl$gene_id,]
stopifnot(rownames(expr.mrna) == ensembl$gene_id)

# Get GO terms for the ensembl genes.
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
attrib = c("ensembl_gene_id", "go_id", "name_1006", "namespace_1003")
filt = c("ensembl_gene_id")

go = getBM(attributes = attrib, filters = filt, values = unique(ensembl$gene_id),
          mart = mart)

# Convert this ensembl ID based list to a GO term based list.
go.cats = unique(go$go_id)
go.domain = split(go, go$namespace_1003)
go.bp = split(go.domain$biological_process$ensembl_gene_id, 
              go.domain$biological_process$go_id)
go.mf = split(go.domain$molecular_function$ensembl_gene_id,
              go.domain$molecular_function$go_id)
go.cc = split(go.domain$cellular_component$ensembl_gene_id,
              go.domain$cellular_component$go_id)
rm(go.domain)

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

get.C.matrix = function(go, genes, min.size, max.size) {

  go.len = sapply(go, length)
  go = go[go.len >= min.size & go.len <= max.size]

  mat = matrix(0, length(genes), length(go), dimnames = list(genes, names(go)))
  for(i in 1:length(go)) {

    mat[match(go[[i]], genes),i] = 1

  } # for(i)
  mat = mat[rowSums(mat) > 0,]

  return(mat)

} # get.C.matrix()

# Replacement safedag function.
safedag = function(object = NULL, ontology = NULL, top = NULL, color.cutoffs = c(0.1,
    0.01, 0.001), filter = 0, max.GOnames = 200)
{
    if (is.null(ontology)) {
        ontology <- Ontology(mget(names(object@global.stat)[1],
            GOTERM)[[1]])
        cat(paste("Ontology unspecified so assumed to be GO.",
            ontology, "\n", sep = ""))
    } else if (!ontology %in% c("GO.CC", "GO.BP", "GO.MF")) {
        stop("ontology must be \"GO.CC\", \"GO.BP\", or \"GO.MF\"",
            call. = FALSE)
	} else {
        ontology = substr(ontology, 4, 5)
    }
	
    if (!is.null(top)) {
        keep2 <- top
        stop <- 0
        names <- names(as.list(get(paste("GO", ontology, "CHILDREN",
            sep = ""))))
        while (stop == 0) {
            parents <- keep2[keep2 %in% names]
            children <- unlist(mget(parents, get(paste("GO",
                ontology, "CHILDREN", sep = ""))))
            stop <- all(children %in% keep2)
            keep2 <- unique(c(keep2, children))
        }
    } else {
        top <- c("GO:0005575", "GO:0008150", "GO:0003674")[c("CC",
            "BP", "MF") %in% ontology]
        keep2 <- NULL
    }
	
    parents <- get(paste("GO", ontology, "PARENTS", sep = ""))
    C.names <- names(object@global.stat)
    keep <- C.names %in% names(as.list(parents))
    cat(paste(sum(keep), " of ", length(C.names), " categories in GO.",
        ontology, "\n", sep = ""))
    if (!is.null(keep2)) {
        keep[!C.names %in% keep2] <- FALSE
        cat(paste(sum(keep), "of", length(C.names), "categories below",
            top, "\n"))
    }
    p <- object@global.pval
    cat(paste("\n", sum(p <= color.cutoffs[3]), "categories with p <=",
        color.cutoffs[3], "\n"))
    if (filter < 3)
        cat(paste("", sum(p <= color.cutoffs[2]), "categories with p <=",
            color.cutoffs[2], "\n"))
    if (filter < 2)
        cat(paste("", sum(p <= color.cutoffs[1]), "categories with p <=",
            color.cutoffs[1], "\n"))
    if (filter)
        keep[p > color.cutoffs[filter]] <- FALSE
    color <- 8 - (p <= color.cutoffs[1]) * 6 + (p <= color.cutoffs[2]) +
        (p <= color.cutoffs[3])
    g <- GOGraph(C.names[keep], parents)
    nda <- makeNodeAttrs(g)
    nodes <- names(nda$fillcolor)
    if (is.null(keep2)) {
        g <- subGraph(nodes[!nodes %in% c("all", "top")], g)
    } else { 
	    g <- subGraph(nodes[nodes %in% keep2], g)
	}
    nda <- makeNodeAttrs(g)
    nodes <- names(nda$fillcolor)
    match <- match(nodes, C.names[keep], nomatch = 0)
    nda$fillcolor[match] <- color[keep]
    nda$fillcolor[!nodes %in% C.names] <- "white"
    if (length(nodes) > max.GOnames) {
        nda$label[nda$label != ""] <- " "
    } else {
        if (is.null(keep2))
            nda$label[nodes == top] <- ontology
    }
    g <- agopen(g, nodeAttrs = nda, name = "whatever")
    par(mfrow = c(1, 1))
    plot(g)

} # safedag()



safe.fxn = function(expr, go.list, y.vec, file.prefix, model = "default",
           Z.mat = NULL) {

  C.mat = getCmatrix(go.list, present.genes = rownames(expr), min.size = 5, 
          max.size = 600)
#  C.mat = get.C.matrix(go = go.list, genes = rownames(expr), min.size = 5, 
#          max.size = 600)

  # Subset the expression again to match the C.matrix.
  expr = expr[C.mat$row.names,]
  print(paste(nrow(expr), "genes &", ncol(C.mat$C.mat.csr), "categories"))
  stopifnot(rownames(expr) == C.mat$row.names)

  results = safe(X.mat = expr, y.vec = y.vec, C.mat = C.mat, local = model,
            alpha = 1.0, error = "FDR.BH", method = "bootstrap", Pi.mat = 10000, parallel = TRUE)

  saveRDS(results, file = paste0(file.prefix, ".rds"))

  pv = results@global.pval
  nm = mget(names(pv), GOTERM, ifnotfound = NA)

  pv = pv[!is.na(nm)]
  nm = nm[!is.na(nm)]

  table = data.frame(GOID = names(pv), Name = sapply(nm, "Term"), pv)
  saveRDS(table, file = paste0(file.prefix, "_table.rds"))

} # safe.fxn()


# Register parallele back end.
registerDoParallel(cores = 19)

# Sex & Diet one-way ANOVA.
expr.res = expr.mrna
y.vec = apply(covar[,2:3], 1, paste, collapse = ":")
safe.fxn(expr = expr.res, go.list = go.bp, y.vec = y.vec,
         file.prefix = "paste0(out_dir, SAFE/DO478_safe_GOBP_sex_diet")
safe.fxn(expr = expr.res, go.list = go.mf, y.vec = y.vec,
         file.prefix = "paste0(out_dir, SAFE/DO478_safe_GOMF_sex_diet")
safe.fxn(expr = expr.res, go.list = go.cc, y.vec = y.vec,
         file.prefix = "paste0(out_dir, SAFE/DO478_safe_GOCC_sex_diet")

# Sex only, regress out diet and gen first.
expr.res = expr.mrna
for(i in 1:nrow(expr.mrna)) {

  mod = lm(expr.mrna[i,] ~ covar[,3:8])
  expr.res[i,] = residuals(mod)

} # for(i)

y.vec = covar[,2]
safe.fxn(expr = expr.res, go.list = go.bp, y.vec = y.vec,
         file.prefix = paste0(out_dir, "DO478_safe_GOBP_sex"))
safe.fxn(expr = expr.res, go.list = go.mf, y.vec = y.vec,
         file.prefix = paste0(out_dir, "DO478_safe_GOMF_sex"))
safe.fxn(expr = expr.res, go.list = go.cc, y.vec = y.vec,
         file.prefix = paste0(out_dir, "DO478_safe_GOCC_sex"))

# Diet only, regress out sex and gen first.
expr.res = expr.mrna
for(i in 1:nrow(expr.mrna)) {

  mod = lm(expr.mrna[i,] ~ covar[,c(2,4:8)])
  expr.res[i,] = residuals(mod)

} # for(i)

y.vec = covar[,3]
safe.fxn(expr = expr.res, go.list = go.bp, y.vec = y.vec,
         file.prefix = paste0(out_dir, "DO478_safe_GOBP_diet"))
safe.fxn(expr = expr.res, go.list = go.mf, y.vec = y.vec,
         file.prefix = paste0(out_dir, "DO478_safe_GOMF_diet"))
safe.fxn(expr = expr.res, go.list = go.cc, y.vec = y.vec,
         file.prefix = paste0(out_dir, "DO478_safe_GOCC_diet"))

# Sex X Diet ANOVA.
# Regress out generation.
expr.res = expr.mrna
for(i in 1:nrow(expr.mrna)) {

  mod = lm(expr.mrna[i,] ~ covar[,4:8])
  expr.res[i,] = residuals(mod)

} # for(i)
y.vec = apply(covar[,2:3], 1, paste, collapse = ":")

safe.fxn(expr = expr.res, go.list = go.bp, y.vec = y.vec,
         file.prefix = paste0(out_dir, "DO478_safe_GOBP_sexXdiet"),
         model = "inter.model")
safe.fxn(expr = expr.res, go.list = go.mf, y.vec = y.vec,
         file.prefix = paste0(out_dir, "DO478_safe_GOMF_sexXdiet"),
         model = "inter.model")
safe.fxn(expr = expr.res, go.list = go.cc, y.vec = y.vec,
         file.prefix = paste0(out_dir, "DO478_safe_GOCC_sexXdiet"),
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
gm = readRDS("safe_GOBP_cat_gene_means.rds")

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




