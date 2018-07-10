expr = expr.res
kegg.list = kegg.paths
y.vec = y.vec
file.prefix = paste0(outdir, "SexAgeInteraction")
model = "inter.model"
Z.mat = NULL
C.mat = getCmatrix(kegg.list, present.genes = rownames(expr), min.size = 5,
        max.size = 600)
expr = expr[C.mat$row.names,]
print(paste(nrow(expr), "genes &", ncol(C.mat$C.mat.csr), "categories"))
stopifnot(rownames(expr) == C.mat$row.names)

#results = safe(X.mat = expr, y.vec = y.vec, C.mat = C.mat, local = model,
#          alpha = 1.0, error = "FDR.BH", method = "bootstrap", Pi.mat = 1000, parallel = TRUE)

X.mat = expr
y.vec = y.vec
C.mat = C.mat
local = model
alpha = 1.0
error = "FDR.BH"
method = "bootstrap"
Pi.mat = 1000
parallel = TRUE

platform = NULL
annotate = NULL
min.size = 2
max.size = Inf
by.gene = FALSE
global = "default"
args.local = NULL
args.global = list(one.sided = FALSE)

epsilon = 10^(-10)
print.it = TRUE


# Custom function fors sex * age interaction.
local.inter.model = function(X.mat, y.vec, ...) {

  sex  = as.numeric(sapply(y.vec, function(x) str_sub(x, str_locate(x, ":")[1]+1,)))
  age = as.numeric(sapply(y.vec, function(x) str_sub(x, ,str_locate(x, ":")[1]-1)))
  local.covar = model.matrix(~ Sex + Age + Sex*Age, data = annot.samples)[,-1]

  return(
    function(data, vector = local.covar, ...) {

      # Sex + Age model
      qry = qr(local.covar[,1:2])
      resid.null = apply(data, 1, function(z) { qr.resid(qry, z) })
      # Sex + Age + Sex*Age model
      qry = qr(local.covar)
      resid.full = apply(data, 1, function(z) { qr.resid(qry, z) })
      return(-ncol(data) * log(colSums(resid.full^2) / colSums(resid.null^2)))
    }
  )
} # local.inter.model()


    df.p <- 0
    C.names <- C.mat$col.names
    C.mat <- C.mat$C.mat.csr

    num.cats <- ncol(C.mat)
    num.genes <- nrow(X.mat)
    num.arrays <- ncol(X.mat)

    if (global == "default")
        global <- "Wilcoxon"
    local.stat <- get(paste("local", local, sep = "."))(X.mat,
        y.vec, args.local)
    u.obs <- local.stat(data = X.mat)
    u.pval <- as.numeric(rep(NA, num.genes))
    if (!is.logical(args.global$one.sided))
        stop("args.global$one.sided is missing or incorrect",
            call. = FALSE)
    global.stat <- get(paste("global", global, sep = "."))(C.mat,
        u.obs, args.global)
    v.obs <- global.stat(u.obs)
    v.pval <- as.numeric(rep(NA, num.cats))


    if (!is.matrix(Pi.mat)) {
        if (Pi.mat == 1) {
            names(u.pval) <- names(u.obs) <- rownames(X.mat)
            names(v.pval) <- names(v.obs) <- C.names
            return(new("SAFE", local = local, local.stat = u.obs,
                local.pval = u.pval, global = global, global.stat = v.obs,
                global.pval = v.pval, error = error, global.error = v.pval,
                alpha = 1, C.mat = C.mat, method = "no resampling"))
        } else {
            num.perms <- Pi.mat
            n = num.arrays

                Pi.mat <- getPImatrix(n = length(y.vec), K = Pi.mat,
                method = method)
                num.perms <- nrow(Pi.mat)
                if (length(y.vec) != ncol(Pi.mat)){
                  stop("Dimensions of Pi.mat and y.vec do not conform",
                  call. = FALSE)
                }
              }
            }
          }

if (method == "bootstrap" | method == "bootstrap.t" | method == "bootstrap.q") {
    if (local %in% c("f.GLM")){
      stop(paste("local = \"", local, "\" can not be used in the bootstrap"),
          call. = FALSE)
    }

    if (global %in% c("Kolmogorov", "Fisher")){
      stop(paste("global = \"", global, "\" cant be used in the bootstrap"),
          call. = FALSE)
    }

    if (error %in% c("FWER.WY", "FDR.YB")){
      stop(paste("error = \"", error, "\" can not be used in the bootstrap"),
          call. = FALSE)
    }

    u.pvalue <- rep(1/num.perms, num.genes)
    u.sum <- u.obs
    u2.sum <- u.obs^2
    null.local <- 0
    emp.p <- rep(1/num.perms, num.cats)
    if (global == "Wilcoxon") {
        C.size <- (rep(1, num.genes) %*% C.mat)[1, ]
        null.global <- (num.genes + 1) * C.size/2
    } else {
      null.global = 0
    }

    v.sum <- v.obs
    v2.sum <- v.obs^2
    if (parallel == FALSE) {
        for (i in 2:num.perms) {
            u <- local.stat(data = X.mat[, Pi.mat[i, ]],
              vector = y.vec[Pi.mat[i, ]], resample = Pi.mat[i,
                ])
            u.sum <- u + u.sum
            u2.sum <- u^2 + u2.sum
            u.pvalue <- u.pvalue + (u * sign(u.obs) <= -epsilon)/num.perms
            v <- global.stat(u)
            v.sum <- v + v.sum
            v2.sum <- v^2 + v2.sum
            emp.p <- emp.p + (v <= (null.global - epsilon))/num.perms
            if (print.it){
              if (trunc(i/100) == i/100){
                cat(paste(i, "bootstrap resamples completed\n"))
              }
            }
        }
    }
    if (method == "bootstrap" | method == "bootstrap.t") {
        emp.p <- 1 - pt(((v.sum/num.perms) - null.global)/sqrt((v2.sum -
            v.sum^2/num.perms)/(num.perms - 1)), df = num.arrays -
            1)
        u.pvalue <- 1 - pt(abs(u.sum/num.perms)/sqrt((u2.sum -
            u.sum^2/num.perms)/(num.perms - 1)), df = num.arrays -
            1)
    }
    if (error == "none") {
        error.p <- as.numeric(rep(NA, num.cats))
        if (is.na(alpha)){
          alpha <- 0.05
        }
    } else {
        error.p <- get(paste("error", error, sep = "."))(t(emp.p))
        if (is.na(alpha)){
            alpha <- 0.1
        }
    }

        names(u.pvalue) <- names(u.obs) <- rownames(X.mat)
        names(error.p) <- names(emp.p) <- names(v.obs) <- C.names
        return(new("SAFE", local = local, local.stat = u.obs,
            local.pval = u.pvalue, global = global, global.stat = v.obs,
            global.pval = emp.p, error = error, alpha = alpha,
            global.error = error.p, C.mat = C.mat, method = method))
