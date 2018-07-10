> safe
function (X.mat, y.vec, C.mat = NULL, Z.mat = NULL, method = "permutation",
    platform = NULL, annotate = NULL, min.size = 2, max.size = Inf,
    by.gene = FALSE, local = "default", global = "default", args.local = NULL,
    args.global = list(one.sided = FALSE), Pi.mat = NULL, error = "FDR.BH",
    parallel = FALSE, alpha = NA, epsilon = 10^(-10), print.it = TRUE,
    ...)
{
    if (!class(X.mat) %in% c("matrix", "data.frame"))
        stop("SAFE v3.* nolonger supports exprSet and other S4 class objects.",
            call. = FALSE)
    if (length(y.vec) != ncol(X.mat))
        stop("Dimensions of X.mat and y.vec do not conform",
            call. = FALSE)
    if (!is.null(Z.mat)) {
        if (length(y.vec) != nrow(Z.mat))
            stop("Dimensions of Z.mat and y.vec do not conform",
                call. = FALSE)
        X.mat <- getXYresiduals(X.mat, y.vec, Z.mat)
        df.p <- X.mat$df.p
        y.vec <- X.mat$y.star
        X.mat <- X.mat$X.star
    }
    else df.p <- 0
    if (is.null(C.mat)) {
        if (is.null(platform) | is.null(annotate)) {
            stop("C.mat or platform & annotate must be specified",
                call. = FALSE)
        }
        else {
            require(platform, character.only = TRUE)
            platform <- sub("[.]db$", "", platform)
            probes <- rownames(X.mat)
            genes <- get(paste(platform, "SYMBOL", sep = ""))
            if (is.null(probes) | (!prod(probes %in% names(as.list(genes))))) {
                stop(paste("Rownames of X.mat do not conform with ",
                  platform, sep = ""), call. = FALSE)
            }
            else genes <- unlist(mget(probes, genes))
            if (substr(annotate[1], 1, 3) == "GO.") {
                cat(paste("Building categories from ", platform,
                  "GO2ALLPROBES\n", sep = ""))
                C.mat <- getCmatrix(keyword.list = as.list(get(paste(platform,
                  "GO2ALLPROBES", sep = ""))), present.genes = probes,
                  GO.ont = substr(annotate, 4, 5), min.size = min.size,
                  max.size = max.size, by.gene = by.gene, gene.names = genes)
                C.names <- C.mat$col.names
                C.mat <- C.mat$C.mat.csr
            }
            else if (annotate == "PFAM") {
                cat(paste("Building categories from ", platform,
                  "PFAM\n", sep = ""))
                C.mat <- getCmatrix(gene.list = as.list(get(paste(platform,
                  "PFAM", sep = ""))), present.genes = probes,
                  min.size = min.size, max.size = max.size, by.gene = by.gene,
                  gene.names = genes)
                C.names <- paste("PFAM:", substr(C.mat$col.names,
                  3, 100), sep = "")
                C.mat <- C.mat$C.mat.csr
            }
            else if (annotate == "REACTOME") {
                require(reactome.db)
                cat(paste("Building categories from reactome.db by ENTREZID"))
                entrez <- as.character(unlist(mget(probes, get(paste(platform,
                  "ENTREZID", sep = "")))))
                entrez[is.na(entrez)] <- paste("NaN", 1:sum(is.na(entrez)))
                C.mat <- getCmatrix(gene.list = as.list(reactomeEXTID2PATHID),
                  present.genes = entrez, min.size = min.size,
                  max.size = max.size, by.gene = by.gene, gene.names = genes)
                C.names <- paste("REACTOME:", C.mat$col.names,
                  sep = "")
                C.mat <- C.mat$C.mat.csr
            }
            else stop(paste("Annotate = \"", annotate, "\" not recognized",
                sep = ""), call. = FALSE)
        }
    }
    else if (class(C.mat) == "matrix") {
        C.names <- colnames(C.mat)
        if (is.null(C.names)) {
            C.names <- paste("Column", 1:ncol(C.mat))
            cat("Warning: Arbitrary columns names generated for categories \n")
        }
        C.mat <- as.matrix.csr(C.mat)
        if (sum(rownames(X.mat) != rownames(C.mat)))
            cat("Warning: gene labels do not match between X.mat and C.mat \n")
    }
    else {
        C.names <- C.mat$col.names
        C.mat <- C.mat$C.mat.csr
    }
    if (nrow(C.mat) != nrow(X.mat)) {
        stop("Dimensions of X.mat and C.mat do not conform",
            call. = FALSE)
    }
    num.cats <- ncol(C.mat)
    num.genes <- nrow(X.mat)
    num.arrays <- ncol(X.mat)
    if (method == "express") {
        if (!local %in% c("default", "Score"))
            cat("WARNING: Only Score statistics can be used in safe-express \n")
        local <- "Score"
        if (!global %in% c("default", "D", "U", "V"))
            cat(paste("WARNING: global = \"", global, "\" not available in safe-express\n",
                sep = ""))
        if (!global %in% c("D", "U", "V"))
            global <- "D"
        if (error %in% c("FWER.WY", "FDR.YB")) {
            cat(paste("WARNING: error = \"", error, "\" not available in safe-express \n",
                sep = ""))
            error <- "none"
        }
        if (!prod(unique(as.matrix(C.mat)) %in% 0:1))
            stop("Soft categories not available in safe-express",
                call. = FALSE)
        if (length(unique(y.vec)) == 2)
            if (!prod(y.vec %in% 0:1)) {
                cat(paste("Warning: y.vec is not (0,1), thus Group 1 ==",
                  y.vec[1], "\n"))
                y.vec <- (y.vec == y.vec[1]) * 1
            }
        C.list <- apply(as.matrix(C.mat) == 1, 2, which)
        if (is.null(args.global$grid.length))
            grid <- 100
        else grid <- args.global$grid.length
        if (is.null(args.global$gamma.thresh))
            gamma <- 0.01
        else gamma <- args.global$gamma.thresh
        if (is.null(args.global$method.bivar))
            mb <- "exact"
        else mb <- args.global$method.bivar
        if (is.null(args.global$swivel))
            sw <- 1e-04
        else sw <- args.global$swivel
        express <- safe.express(X.mat = X.mat, y.vec = y.vec,
            C.list = C.list, df.penalty = df.p, method = global,
            gamma.thresh = gamma, method.bivar = mb, grid = grid,
            swivel = sw, parallel = parallel)
        if (global == "D") {
            global.pval <- unlist(express$global.D[, 2])
            global.stat <- unlist(express$global.D[, 1])
        }
        if (global == "V") {
            global.pval <- unlist(express$global.V[, 2])
            global.stat <- unlist(express$global.V[, 1])
        }
        if (global == "U") {
            global.pval <- unlist(express$global.U[, 2])
            global.stat <- unlist(express$global.U[, 1])
        }
        if (error == "none") {
            global.error <- as.numeric(rep(NA, num.cats))
            if (is.na(alpha))
                alpha <- 0.05
        }
        else {
            global.error <- get(paste("error", error, sep = "."))(t(global.pval))
            if (is.na(alpha))
                alpha <- 0.1
        }
        local.stat <- express$local.stat
        local.pval <- pchisq(local.stat^2, 1, lower.tail = F)
        names(local.stat) <- names(local.pval) <- rownames(X.mat)
        names(global.stat) <- names(global.pval) <- names(global.error) <- C.names
        return(new("SAFE", local = local, local.stat = local.stat,
            local.pval = local.pval, global = global, global.stat = global.stat,
            global.pval = global.pval, error = error, alpha = alpha,
            global.error = global.error, C.mat = C.mat, method = method))
    }
    if (local == "default") {
        if (length(unique(y.vec)) == 2) {
            local <- "t.Student"
        }
        else if (class(y.vec) == "character") {
            local <- "f.ANOVA"
        }
        else local <- "t.LM"
    }
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
    if (is.null(Pi.mat))
        if (substr(method, 1, 4) == "boot")
            Pi.mat <- 200
        else Pi.mat <- 1000
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
            if (local %in% c("t.paired")) {
                if (method != "permutation")
                  count <- factorial(n/2)
                else count <- 2^(n/2 - 1)
                if (count < num.perms)
                  cat(paste("Warning: only", round(count), "unique resamples exist\n"))
                Pi.mat <- getPImatrix(block.vec = y.vec, K = num.perms,
                  method = method)
            }
            else if (local %in% c("t.Student", "t.Welch") & method ==
                "permutation") {
                count <- choose(n, table(y.vec)[1])
                if (count < num.perms) {
                  cat(paste("Warning: only", round(count), "unique resamples exist\n",
                    "         switching to exhaustive permutation\n"))
                  Pi.mat <- getPIcomplete(y.vec)
                  num.perms <- nrow(Pi.mat)
                }
                else Pi.mat <- getPImatrix(y.vec = y.vec, K = num.perms,
                  method = method)
            }
            else if (local %in% c("f.ANOVA", "t.Student", "t.Welch")) {
                n2 <- table(y.vec)
                n3 <- 1
                for (i in 1:length(n2)) n3 <- n3 * sum(choose(n2[i],
                  k = 2:n2[i]) * choose(n2[i] - 1, k = 1:(n2[i] -
                  1)))
                if (method != "permutation")
                  count <- n3
                else count <- exp(lgamma(n + 1) - sum(lgamma(n2 +
                  1)))
                if (count < num.perms)
                  cat(paste("Warning: only", round(count), "unique resamples exist\n"))
                Pi.mat <- getPImatrix(y.vec = y.vec, K = Pi.mat,
                  method = method)
            }
            else if (local %in% c("t.LM", "z.COXPH")) {
                if (method != "permutation")
                  count <- sum(choose(n, k = 2:n) * choose(n -
                    1, k = 1:(n - 1)))
                else count <- factorial(n)
                if (count < num.perms)
                  cat(paste("Warning: only", round(count), "unique resamples exist\n"))
                Pi.mat <- getPImatrix(n = length(y.vec), K = Pi.mat,
                  method = method)
            }
            else {
                Pi.mat <- getPImatrix(n = length(y.vec), K = Pi.mat,
                  method = method)
            }
        }
    }
    else {
        num.perms <- nrow(Pi.mat)
        if (length(y.vec) != ncol(Pi.mat))
            stop("Dimensions of Pi.mat and y.vec do not conform",
                call. = FALSE)
    }
    if (method == "permutation") {
        u.pvalue <- rep(1/num.perms, num.genes)
        if (error %in% c("none", "FWER.Bonf", "FDR.BH")) {
            error.p <- v.pval
            emp.p <- rep(1/num.perms, num.cats)
            if (parallel == FALSE) {
                for (i in 2:num.perms) {
                  u <- local.stat(data = X.mat[, Pi.mat[i, ]],
                    resample = Pi.mat[i, ])
                  u.pvalue <- u.pvalue + (abs(u) >= (abs(u.obs) +
                    epsilon))/num.perms
                  v <- global.stat(u)
                  emp.p <- emp.p + (v >= (v.obs + epsilon))/num.perms
                  if (print.it)
                    if (trunc(i/100) == i/100)
                      cat(paste(i, "permutations completed\n"))
                }
            }
            else {
                require(foreach)
                require(doRNG)
                if (print.it)
                  cat(paste("Permutations split across", getDoParWorkers(),
                    "cores\n"))
                parallel.p <- foreach(i = 2:num.perms, .combine = "+",
                  .inorder = FALSE) %dorng% {
                  u <- local.stat(data = X.mat[, Pi.mat[i, ]],
                    resample = Pi.mat[i, ])
                  u.frac <- (abs(u) >= (abs(u.obs) + epsilon))
                  v <- global.stat(u)
                  emp.frac <- (v >= (v.obs + epsilon))
                  c(u.frac, emp.frac)
                }
                u.pvalue <- u.pvalue + parallel.p[1:length(u.obs)]/num.perms
                emp.p <- emp.p + parallel.p[(length(u.obs) +
                  1):(length(u.obs) + length(v.obs))]/num.perms
            }
            if (error != "none") {
                error.p <- get(paste("error", error, sep = "."))(t(emp.p))
            }
            if (is.na(alpha))
                alpha <- 0.05
        }
        else {
            if (parallel == FALSE) {
                V.mat <- matrix(0, num.perms, num.cats)
                V.mat[1, ] <- v.obs
                for (i in 2:num.perms) {
                  u <- local.stat(data = X.mat[, Pi.mat[i, ]],
                    resample = Pi.mat[i, ])
                  u.pvalue <- u.pvalue + (abs(u) >= (abs(u.obs) +
                    epsilon))/num.perms
                  V.mat[i, ] <- global.stat(u)
                  if (print.it)
                    if (trunc(i/100) == i/100)
                      cat(paste(i, "permutations completed\n"))
                }
            }
            else stop(paste("error = \"", error, "\" can not be used with parallel = TRUE"),
                call. = FALSE)
            P.mat <- (num.perms + 1 - apply(V.mat, 2, rank, ties.method = "min"))/num.perms
            error.p <- get(paste("error", error, sep = "."))(P.mat)
            emp.p <- P.mat[1, , drop = TRUE]
            if (is.na(alpha))
                alpha <- 0.1
        }
        names(u.pvalue) <- names(u.obs) <- rownames(X.mat)
        names(error.p) <- names(emp.p) <- names(v.obs) <- C.names
        return(new("SAFE", local = local, local.stat = u.obs,
            local.pval = u.pvalue, global = global, global.stat = v.obs,
            global.pval = emp.p, error = error, alpha = alpha,
            global.error = error.p, C.mat = C.mat, method = method))
    }
    else if (method == "bootstrap" | method == "bootstrap.t" |
        method == "bootstrap.q") {
        if (local %in% c("f.GLM"))
            stop(paste("local = \"", local, "\" can not be used in the bootstrap"),
                call. = FALSE)
        if (global %in% c("Kolmogorov", "Fisher"))
            stop(paste("global = \"", global, "\" cant be used in the bootstrap"),
                call. = FALSE)
        if (error %in% c("FWER.WY", "FDR.YB"))
            stop(paste("error = \"", error, "\" can not be used in the bootstrap"),
                call. = FALSE)
        u.pvalue <- rep(1/num.perms, num.genes)
        u.sum <- u.obs
        u2.sum <- u.obs^2
        null.local <- 0
        emp.p <- rep(1/num.perms, num.cats)
        if (global == "Wilcoxon") {
            C.size <- (rep(1, num.genes) %*% C.mat)[1, ]
            null.global <- (num.genes + 1) * C.size/2
        }
        else null.global = 0
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
                if (print.it)
                  if (trunc(i/100) == i/100)
                    cat(paste(i, "bootstrap resamples completed\n"))
            }
        }
        else {
            require(foreach)
            require(doRNG)
            if (print.it)
                cat(paste("Bootstrap resamples split across",
                  getDoParWorkers(), "cores\n"))
            parallel.p <- foreach(i = 2:num.perms, .combine = "+",
                .inorder = FALSE) %dorng% {
                u <- local.stat(data = X.mat[, Pi.mat[i, ]],
                  vector = y.vec[Pi.mat[i, ]], resample = Pi.mat[i,
                    ])
                u.frac <- (u * sign(u.obs) <= -epsilon)
                v <- global.stat(u)
                emp.frac <- (v <= (null.global - epsilon))
                c(u, u^2, u.frac, v, v^2, emp.frac)
            }
            split.u <- length(u.obs)
            split.v <- length(v.obs)
            u.sum <- u.sum + parallel.p[1:split.u]
            u2.sum <- u2.sum + parallel.p[(split.u + 1):(2 *
                split.u)]
            u.pvalue <- u.pvalue + parallel.p[(2 * split.u +
                1):(3 * split.u)]/num.perms
            v.sum <- v.sum + parallel.p[(3 * split.u + 1):(3 *
                split.u + split.v)]
            v2.sum <- v2.sum + parallel.p[(3 * split.u + split.v +
                1):(3 * split.u + 2 * split.v)]
            emp.p <- emp.p + parallel.p[(3 * split.u + 2 * split.v +
                1):(3 * split.u + 3 * split.v)]/num.perms
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
            if (is.na(alpha))
                alpha <- 0.05
        }
        else {
            error.p <- get(paste("error", error, sep = "."))(t(emp.p))
            if (is.na(alpha))
                alpha <- 0.1
        }
        names(u.pvalue) <- names(u.obs) <- rownames(X.mat)
        names(error.p) <- names(emp.p) <- names(v.obs) <- C.names
        return(new("SAFE", local = local, local.stat = u.obs,
            local.pval = as.numeric(u.pvalue), global = global,
            global.stat = v.obs, global.pval = emp.p, error = error,
            alpha = alpha, global.error = error.p, C.mat = C.mat,
            method = method))
    }
    else stop(paste("method = \"", method, "\" is not recognized",
        sep = ""), call. = FALSE)
}
