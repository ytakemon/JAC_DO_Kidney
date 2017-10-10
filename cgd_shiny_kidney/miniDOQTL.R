add.sig.thr = function(sig.thr, sig.col = "red", chrsum) {

  # If sig.thr is a matrix, then we have separate autosome and X thresholds.
  if(is.matrix(sig.thr)) {

    # Make sure that the thresholds and colors have the same length.
    if(length(sig.col) < nrow(sig.thr)) {
      stop(paste("The number of colors is less than the number of thresholds.",
           "Please provide as many colors as there are rows in the sig.thr matrix."))
    } # if(length(sig.col) < nrow(sig.thr))

    # Get the coordinates to the end of the autosomes.
    old.warn = options("warn")$warn
    options(warn = -1)    
    tmp = as.numeric(names(chrsum))
    auto.len = max(chrsum[!is.na(tmp)])
    X.len = max(chrsum[is.na(tmp)])
    options(warn = old.warn)

    for(i in 1:nrow(sig.thr)) {
      lines(x = c(0, auto.len), y = rep(sig.thr[i,"A"], 2), lwd = 2,
            col = sig.col[i])
      lines(x = c(auto.len, X.len), y = rep(sig.thr[i,"X"], 2), lwd = 2,
             col = sig.col[i])
    } # for(i)

  } else {

    # Make sure that the thresholds and colors have the same length.
    if(length(sig.col) < length(sig.thr)) {
      stop(paste("The number of colors is less than the number of thresholds.",
           "Please provide as many colors as there are thresholds in sig.thr."))
    } # if(length(sig.col) < nrow(sig.thr))

    abline(h = sig.thr, lwd = 2, col = sig.col)

  } # else

} # add.sig.thr


plot.doqtl = function(x, stat.name = c("lod", "neg.log10.p"),  sig.thr = NULL, 
                      sig.col = "red", ...) {
  
  doqtl = x
  chrlen = structure(c(195.471971, 182.113224, 160.03968, 156.508116, 151.834684,
                       149.736546, 145.441459, 129.401213, 124.59511, 130.694993, 122.082543,
                       120.129022, 120.421639, 124.902244, 104.043685, 98.207768, 94.987271,
                       90.702639, 61.431566, 171.031299, 91.744698, 0.016299), 
                       .Names = c(as.character(1:19), "X", "Y", "M"))
  stat.name = match.arg(stat.name)

  # Get the call and arguments.
  call = match.call()
  lod = doqtl$lod$A
  if(any(names(doqtl$lod) == "X")) {
    lod = rbind(doqtl$lod$A, doqtl$lod$X)
  } # if(length(lod) > 1)

  # Get chr lengths and locations.  Create an x-axis based on Genome Mb.
  if(max(lod[,3], na.rm = TRUE) > 200) {
    lod[,3] = lod[,3] * 1e-6
  } # if(max(lod[,3]) > 200)

  mb = lod[,3]
  gmb = mb
  unique.chr = as.character(unique(lod[,2]))
  old.warn = options("warn")$warn
  options(warn = -1)
  unique.chr = unique.chr[order(as.numeric(unique.chr))]
  options(warn = old.warn)
  chrlen = chrlen[names(chrlen) %in% unique.chr]

  if(max(chrlen) > 200) {
    chrlen = chrlen * 1e-6
  } # if(max(chrlen) > 200)

  # Get the cumulative sum of the Chr lengths.
  chrsum = cumsum(chrlen)
  chrsum = c(0, chrsum)
  chrmid = chrsum[-length(chrsum)] + (diff(chrsum) * 0.5)

  # Add the preceeding chromosome lengths to each SNP position.
  for(c in 2:length(unique.chr)) {
    rows = which(lod[,2] == unique.chr[c])
    gmb[rows] = gmb[rows] + chrsum[c]
  } # for(c)

  # Make the basic plot.
  plot.column = which(colnames(lod) == stat.name)
  if(length(plot.column) == 0) {
    stop(paste("The stat.name of", stat.name, "was not found in the column names",
         "of the DOQTL object. Please verify that stat.name contains one of the",
         "column names in the DOQTL object."))
  } # if(length(plot.column) == 0)

  par(font = 2, font.lab = 2, font.axis = 2, las = 1, xaxs = "i",
      plt = c(0.12, 0.95, 0.05, 0.88))

  if("ylim" %in% names(call)) {
    plot(gmb, lod[,plot.column], col = 0, xlab = "", xaxt = "n", ylab = stat.name, ...)
  } else {
    if(!missing(sig.thr)) {
      plot(gmb, lod[,plot.column], col = 0, xlab = "", ylab = stat.name,
           xaxt = "n", ylim = c(0, max(lod[,plot.column], sig.thr, na.rm = TRUE) * 1.05), ...)
    } else {
      plot(gmb, lod[,plot.column], col = 0, xlab = "", ylab = stat.name,
           xaxt = "n", ylim = c(0, max(lod[,plot.column], na.rm = TRUE) * 1.05), ...)
    } # else
  } # else

  lod = cbind(lod, gmb)
  lod = split(lod, lod[,2])
  usr = par("usr")
  rect(chrsum[2 * 1:(length(chrsum) * 0.5) - 1], usr[3],
       chrsum[2 * 1:(length(chrsum) * 0.5)], usr[4], col = rgb(0.8,0.8,0.8),
       border = NA)
  rect(usr[1], usr[3], usr[2], usr[4], border = 1)
  text(chrmid, 0.95 * usr[4], names(chrsum)[-1])
  lapply(lod, function(z) { points(z$gmb, z[,plot.column], type = "l", lwd = 2)})

  if(!is.null(sig.thr)) {
    add.sig.thr(sig.thr = sig.thr, sig.col = sig.col, chrsum = chrsum)
  } # if(!is.null(sig.thr))

} # plot.doqtl() 

bayesint = function(qtl, chr, prob = 0.95, expandtomarkers = FALSE) {
  
  if(missing(qtl)) {
    stop("bayesint: The qtl cannot be null. Please supply a QTL object.")
  } # if(is.null(qtl))
  
  if(missing(chr)) {
    stop(paste("bayesint: The chromosome cannot be null."))
  } else if(!chr %in% c(1:19, "X")) {
    stop(paste("bayesint: The chromosome must be 1 to 19 or X."))
  } # else if
  
  if(prob < 0 | prob > 1) {
    stop(paste("bayesint: The probability must between 0 and 1."))
  } # if(prob < 0 | prob > 1)
  
  old.warn = options("warn")$warn
  options(warn = -1)
  if(!is.na(as.numeric(chr))) {
    qtl = qtl$lod$A
  } else {
    qtl = qtl$lod$X
  } # else
  options(warn = old.warn)
  
  qtl[,1] = as.character(qtl[,1])
  qtl[,2] = as.character(qtl[,2])
  qtl[,3] = as.numeric(qtl[,3])
  qtl[,7] = as.numeric(qtl[,7])
  qtl = qtl[qtl[,2] == chr,]
  pos = qtl[,3]
  
  if(any(is.na(pos))) {
    remove = which(is.na(pos))
    qtl = qtl[-remove,]
    pos = pos[-remove]
  } # if(any(is.na(pos)))
  
  # Make a set of 10,000 intervals so that we can integrate numerically.
  breaks = approx(x = pos, y = 10^qtl[,7], xout = seq(pos[1], pos[length(pos)],
                                                      length.out = 1e5))
  widths  = diff(breaks$x)
  heights = breaks$y[-1] + breaks$y[-length(breaks$y)]
  trapezoids = 0.5 * heights * widths
  # Normalize the area to 1.
  trapezoids = trapezoids / sum(trapezoids)
  # This code copied from the R/qtl bayesint() function by Karl Broman.
  ord = order(breaks$y[-length(breaks$y)], decreasing = TRUE)
  wh  = min(which(cumsum(trapezoids[ord]) >= prob))
  int = range(ord[1:wh])
  # Find the left & right SNP.
  left.snp  = c(NA, qtl[1,2], breaks$x[int][1], 
                approx(qtl[,3], qtl[,4], breaks$x[int][1])$y,
                approx(qtl[,3], qtl[,5], breaks$x[int][1])$y,
                approx(qtl[,3], qtl[,6], breaks$x[int][1])$y,
                approx(qtl[,3], qtl[,7], breaks$x[int][1])$y)
  max.snp   = qtl[which.max(qtl[,7]),]
  right.snp = c(NA, qtl[1,2], breaks$x[int][2], 
                approx(qtl[,3], qtl[,4], breaks$x[int][2])$y,
                approx(qtl[,3], qtl[,5], breaks$x[int][2])$y,
                approx(qtl[,3], qtl[,6], breaks$x[int][2])$y,
                approx(qtl[,3], qtl[,7], breaks$x[int][2])$y)
  if(expandtomarkers) {
    # Find the left & right SNP.
    left.snp  = qtl[max(which(breaks$x[int][1] >= qtl[,3])),]
    max.snp   = qtl[which.max(qtl[,7]),]
    right.snp = qtl[min(which(breaks$x[int][2] <= qtl[,3])),] 
  } # if(expandtomarkers)
  retval = rbind(left.snp, max.snp, right.snp)
  retval[,3] = round(as.numeric(retval[,3]), digits = 6)
  retval[,4] = round(as.numeric(retval[,4]), digits = 6)
  retval[,5] = round(as.numeric(retval[,5]), digits = 6)
  retval[,6] = round(as.numeric(retval[,6]), digits = 6)
  retval$lod = as.numeric(retval[,7])
  return(retval)
} # bayesint()

coefplot = function(doqtl, chr = 1, stat.name = "LOD", conf.int = TRUE, legend = TRUE,
                    colors = "DO", sex, ...) {
  
  old.par = par(no.readonly = TRUE)
  
  cross = attr(doqtl, "cross")
  if(is.null(cross)) {
    if(colors[1] == "DO") {    
      colors = structure(list(CC_Designation = c("A", "B", "C", "D", "E", "F",
                                                 "G", "H"), Strain = c("A/J", "C57BL/6J", "129S1/SvImJ", "NOD/ShiLtJ",
                                                                       "NZO/HlLtJ", "CAST/EiJ", "PWK/PhJ", "WSB/EiJ"), R_Color = c("#F0F000",
                                                                                                                                   "#808080", "#F08080", "#1010F0", "#00A0F0", "#00A000", "#F00000",
                                                                                                                                   "#9000E0")), .Names = c("CC_Designation", "Strain", "R_Color"
                                                                                                                                   ), row.names = c(NA, -8L), class = "data.frame")
    } else if(colors[1] == "HS") {
      colors = hs.colors
    } # else if(colors[1] == "HS")
  } else {
    if(cross == "DO") {    
      colors = do.colors
    } else if(cross == "HS") {
      colors = hs.colors
    } # else if(cross == "HS")
  } # else
  
  num.founders = nrow(colors)
  call = match.call()
  
  # Keep only the founder coefficients from the coef.matrix.
  lod = NULL
  coef = NULL
  if(chr == "X") {
    if(missing(sex)) {
      stop("Sex (either M or F) must be specified on X chromosome.")
    } # if(missing(sex))
    lod  = doqtl$lod$X
    coef = doqtl$coef$X
    if(sex == "F") {
      columns = match(paste("F", colors[,1], sep = "."), colnames(coef))
    } else {
      columns = match(paste("M", colors[,1], sep = "."), colnames(coef))
    } # else
    columns = columns[!is.na(columns)]
    coef = coef[,c(1,columns)]
    colnames(coef)[1] = "A"
    colnames(coef) = sub("^[MF]\\.", "", colnames(coef))
    coef = coef[rownames(coef) %in% lod[,1],]
    # Center the coefficient values.
    coef[,2:ncol(coef)] = coef[,2:ncol(coef)] + coef[,1]
    coef = coef - rowMeans(coef)
  } else {
    lod = doqtl$lod$A
    lod = lod[lod[,2] == chr,]
    intercept = doqtl$coef$A[,1]
    coef = doqtl$coef$A[,(ncol(doqtl$coef$A)-num.founders+1):ncol(doqtl$coef$A)]
    coef[,1] = intercept
    colnames(coef)[1] = "A"
    coef = coef[rownames(coef) %in% lod[,1],]
    # Center the coefficient values.
    coef[,2:ncol(coef)] = coef[,2:ncol(coef)] + coef[,1]
    coef = coef - rowMeans(coef)
  } # else 
  # Verify that the SNP IDs in the lod & coef matrices match.
  if(!all(lod[,1] == rownames(coef))) {
    stop(paste("The SNP IDs in column 1 of the qtl data frame must match",
               "the SNP IDs in the rownames of the coef matrix."))
  } # if(!all(lod[,1] == rownames(coef)))
  # Verify that the coefficient column names are in the colors.
  if(!all(colnames(coef) %in% colors[,1])) {
    stop(paste("The founder names in the colnames of the coefficient matrix",
               "must be in column 1 of the colors matrix."))
  } # if(!all(colnames(coef) %in% colors[,1]))
  # Convert the chromosome locations to Mb.
  if(max(lod[,3], na.rm = TRUE) > 200) {
    lod[,3] = lod[,3] * 1e-6
  } # if(max(lod[,3], na.rm = TRUE) > 200)
  # Set the layout to plot the coefficients on top and the p-values on the 
  # bottom.
  layout(mat = matrix(1:2, 2, 1), heights = c(0.66666667, 0.3333333))
  par(font = 2, font.lab = 2, font.axis = 2, las = 1, plt =
        c(0.12, 0.99, 0, 0.85), xaxs = "i", lwd = 2)
  # Plot the coefficients.
  plot(lod[,3], coef[,colors[1,1]], type = "l", col = colors[1,3], lwd = 2,
       ylim = c(min(coef, na.rm = TRUE), max(coef * 2, na.rm = TRUE)), xlab = 
         paste("Chr", chr), ylab = "Founder Effects", axes = FALSE, ...)
  abline(v = 0:20 * 10, col = "grey80")
  for(i in 1:nrow(colors)) {
    points(lod[,3], coef[,colors[i,1]], type = "l", col = colors[i,3],
           lwd = 2)
  } # for(i)
  # Draw a legend for the founder names and colors.
  if(legend) {
    legend.side = "topleft"
    if(which.max(lod[,7]) < nrow(lod) / 2) {
      legend.side = "topright"
    } # if(which.max(apply(coef, 1, max)) < nrow(lod) / 2)
    legend(legend.side, colors[,2], col = colors[,3], lty = 1, lwd = 2,
           x.intersp = 0.75, y.intersp = 0.75, bg = "white", cex = 0.8)
  } # if(legend)
  # Add the axis.
  axis(2)
  # Plot a rectangle around the plot.
  par(xpd = NA)
  usr = par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], lwd = 2)
  par(xpd = FALSE)
  # Plot the mapping statistic.
  par(plt = c(0.12, 0.99, 0.35, 1))
  # Single chromosome plot.
  plot(lod[,3], lod[,7], type = "l", lwd = 2, xlab = "",
       ylab = stat.name, ...)
  abline(v = 0:20 * 10, col = "grey80")
  points(lod[,3], lod[,7], type = "l", lwd = 2)
  # Shade the confidence interval.
  if(conf.int) {
    interval = bayesint(doqtl, chr = chr)
    usr = par("usr")
    rect(interval[1,3], usr[3], interval[3,3], usr[4], col = rgb(0,0,1,0.1), 
         border = NA)
  } # if(!is.na(conf.int))
  mtext(paste("Chr", chr), 1, 2)
  usr = par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], lwd = 2)
  par(old.par)
} # coefplot()
