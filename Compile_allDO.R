# Find, clean, and combine Aging DO cross-sectional and longitudinal data
# Note: Excluding HDO data

library(rhdf5)
library(dplyr)
library(stringr)

# Get JAC data
# Data: /hpcdata/gac/derived/CGD_DO_Genoprobs/MegaMUGA_hap_probs_v6.h5
# Readme: /hpcdata/gac/derived/CGD_DO_Genoprobs/MegaMUGA_hap_probs_v6.README.txt
# Samples: /hpcdata/gac/raw/JAC_DO_Phenotypes/mouse_info/DO aging study mouse info.xlsx
h5 <- "/hpcdata/gac/derived/CGD_DO_Genoprobs/MegaMUGA_hap_probs_v6.h5"
h5ls(h5)
jac <- h5read(file = h5, name = "/JAC") #probs are naked, need to populate dimnames
dimnames(jac$probs) <- list(jac$samples, jac$founders, jac$markers)
probs <- jac$probs # 1536 x 8 x 68268

# get sample annotation
samples <- read.delim("/hpcdata/gac/raw/JAC_DO_Phenotypes/mouse_info/DO aging study mouse info.txt")
rownames(samples) <- samples$Mouse.ID
# need generation info
key <- read.delim("/hpcdata/gac/raw/JAC_DO_Phenotypes/mouse_info/archive/mouse.key.old.txt")
rownames(key) <- key$Mouse.ID
stopifnot(identical(rownames(key), rownames(samples)))
samples$Generation <- key$Generation

# fix probs naming and replce with proper format
fix <- as.data.frame(rownames(probs))
fix$a <- str_sub(fix[,1],4)
fix$b <- str_split_fixed(fix$a, "[.]", n = 2)[,1]
fix$c <- str_split_fixed(fix$a, "[.]", n = 2)[,2]
fix$d <- str_extract(fix$c, "[[:digit:]]+")
fix$e <- str_pad(fix$d, width = 4, side = "left", pad = 0)
fix$f <- str_c(fix$b, fix$e, sep="-")
fix$dup <- (duplicated(fix$f) | duplicated(fix$f, fromLast = TRUE))
fix <- fix[which((fix$dup %in% FALSE) & (fix$b %in% "DO")) , ] #nrow = 1123
sub <- probs[rownames(probs) %in% fix[,1],,] # dim = 1123 x 8 x 68268
stopifnot(any(rownames(sub) == fix[,1]))
rownames(sub) <- fix$f

# Get genoprobs - sample intersection
samples <- samples[rownames(samples) %in% rownames(sub),] # dim 1123 x 17
# drop unused level ie. Harrison
samples$Generation <- droplevels(samples$Generation)

# Get all urine chemistry
Ucross <- read.delim("/hpcdata/gac/raw/JAC_DO_Phenotypes/urine_chemistry/urine.chemistry.cross.txt", stringsAsFactors = FALSE)
Ulong <- read.delim("/hpcdata/gac/raw/JAC_DO_Phenotypes/urine_chemistry/urine.chemistry.long.txt", stringsAsFactors = FALSE)

# Clean samples and names
Ucross$dup <- (duplicated(Ucross$Mouse.ID) | duplicated(Ucross$Mouse.ID, fromLast = TRUE))
Ucross <- Ucross[Ucross$dup == FALSE, ]
rownames(Ucross) <- Ucross$Mouse.ID
Ucross <- Ucross[,c("Mouse.ID", "Age.Urine.Chem", "cr.u", "mg.u", "mg.cr.u",
                    "ma.u", "ma.cr.u", "phs.u", "phs.cr.u")]

#Ucross data is alrady numeric, but needs to be forced
Ulong$name <- str_sub(Ulong$g, 1, 7)
Ulong$dup <- (duplicated(Ulong$name) | duplicated(Ulong$name, fromLast = TRUE))
Ulong <- Ulong[Ulong$dup == FALSE, ]
rownames(Ulong) <- Ulong$name
Ulong$Mouse.ID <- rownames(Ulong)
Ulong <- Ulong[, c(50,11:45)]
Ulong[,2:36] <- sapply(Ulong[, 2:36], as.numeric)



> str(Ulong)
'data.frame':	593 obs. of  36 variables:
 $ Mouse.ID                       : chr  "DO-0481" "DO-0482" "DO-0483" "DO-0484" ...
 $ X28.wk.Cr.urine..mg.dL.        : chr  "66.03" "95.42" "103.17" "124.86" ...
 $ X28.wk.Mg.urine..mg.dL.        : chr  "25.02" "66.32" "42.3" "89.97" ...
 $ X28.wk.Mg.urine..mg.grcreat.   : chr  "378.92" "695.03" "410.00" "720.57" ...
 $ X28.wk.MA.urine..mg.dL.        : chr  "0" "BDR" "BDR" "0.6" ...
 $ X28.wk.MA.urine..mg.grcreat.   : chr  "0" "BDR" "BDR" "4.81" ...
 $ X28.wk.Phs..mg.dL.             : chr  "117" "240" "360" "590" ...
 $ X28.wk.Phs..mg.grcreat.        : num  1772 2515 3489 4725 3478 ...
 $ X54.wk.Cr.urine..mg.dL.        : num  52 71.3 61.8 NA 117.8 ...
 $ X54.wk.Mg.urine..mg.dL.        : chr  "40.27" "72.49" "46.86" NA ...
 $ X54.wk.Mg.urine..mg.grcreat.   : chr  "775.17" "1017.12" "757.88" NA ...
 $ X54.wk.MA.urine..mg.dL.        : chr  "BDR" "BDR" "BDR" NA ...
 $ X54.wk.MA.urine..mg.grcreat.   : chr  "BDR" "BDR" "BDR" NA ...
 $ X54.wk.Phs..mg.dL.             : num  47.5 194.7 312.4 NA 281.7 ...
 $ X54.wk.Phs..mg.grcreat.        : num  914 2732 5053 NA 2390 ...
 $ X80.week.Cr.urine..mg.dL.      : chr  NA "52.77" "85.03" NA ...
 $ X80.week.Mg.urine..mg.dL.      : num  NA 79.5 83.8 NA 130.3 ...
 $ X80.week.Mg.urine..mg.grcreat. : chr  NA "1505.78" "985.89" NA ...
 $ X80.week.MA.urine..mg.dL.      : chr  NA "BDR" "26.7" NA ...
 $ X80.week.MA.urine..mg.grcreat. : chr  NA "BDR" "314.0" NA ...
 $ X80.week.Phs..mg.dL.           : num  NA 231 240 NA 244 ...
 $ X80.weelkPhs..mg.grcreat.      : chr  NA "4377.49" "2819.01" NA ...
 $ X106..week.Cr.urine..mg.dL.    : chr  NA "30.78" "66.91" NA ...
 $ X106.week.Mg.urine..mg.dL.     : num  NA 56.4 78.4 NA 140 ...
 $ X106.week.Mg.urine..mg.grcreat.: chr  NA "1832.4" "1171.7" NA ...
 $ X106.week.MA.urine..mg.dL.     : chr  NA "1.90" "18.30" NA ...
 $ X106.week.MA.urine..mg.grcreat.: chr  NA "61.7" "273.5" NA ...
 $ X106.week.Phs..mg.dL.          : num  NA 53.7 178.5 NA 211.1 ...
 $ X106.weelkPhs..mg.grcreat.     : chr  NA "1744.6" "2667.8" NA ...
 $ X134..week.Cr.urine..mg.dL.    : num  NA NA NA NA NA ...
 $ X134.week.Mg.urine..mg.dL.     : num  NA NA NA NA NA ...
 $ X134.week.Mg.urine..mg.grcreat.: num  NA NA NA NA NA ...
 $ X134.week.MA.urine..mg.dL.     : chr  NA NA NA NA ...
 $ X134.week.MA.urine..mg.grcreat.: chr  NA NA NA NA ...
 $ X134.week.Phs..mg.dL.          : num  NA NA NA NA NA ...
 $ X134.weelkPhs..mg.grcreat.     : num  NA NA NA NA NA ...





# Prepare data for rbinding
Ucross_prep <- Ucross[c("Mouse.ID")]
Ucross_prep$cr.u.6 <- NA
Ucross_prep$cr.u.12 <- NA
Ucross_prep$cr.u.18 <- NA
Ucross_prep$mg.u.6 <- NA
Ucross_prep$mg.u.12 <- NA
Ucross_prep$mg.u.18 <- NA
Ucross_prep$mg.cr.u.6 <- NA
Ucross_prep$mg.cr.u.12 <- NA
Ucross_prep$mg.cr.u.18 <- NA
Ucross_prep$ma.u.6 <- NA
Ucross_prep$ma.u.12 <- NA
Ucross_prep$ma.u.18 <- NA
Ucross_prep$ma.cr.u.6 <- NA
Ucross_prep$ma.cr.u.12 <- NA
Ucross_prep$ma.cr.u.18 <- NA
Ucross_prep$phs.u.6 <- NA
Ucross_prep$phs.u.12 <- NA
Ucross_prep$phs.u.18 <- NA
Ucross_prep$phs.cr.u.6 <- NA
Ucross_prep$phs.cr.u.12 <- NA
Ucross_prep$phs.cr.u.18 <- NA
Ucross6 <- Ucross[Ucross$Age.Urine.Chem == 6,]
Ucross12 <- Ucross[Ucross$Age.Urine.Chem == 12,]
Ucross18 <- Ucross[Ucross$Age.Urine.Chem == 18,]

Ucross_prep$cr.u.6[which(rownames(Ucross_prep) %in% rownames(Ucross6))] <- Ucross6$cr.u
Ucross_prep$cr.u.12[which(rownames(Ucross_prep) %in% rownames(Ucross12))] <- Ucross12$cr.u
Ucross_prep$cr.u.18[which(rownames(Ucross_prep) %in% rownames(Ucross18))] <- Ucross18$cr.u
Ucross_prep$mg.u.6[which(rownames(Ucross_prep) %in% rownames(Ucross6))] <- Ucross6$mg.u
Ucross_prep$mg.u.12[which(rownames(Ucross_prep) %in% rownames(Ucross12))] <- Ucross12$mg.u
Ucross_prep$mg.u.18[which(rownames(Ucross_prep) %in% rownames(Ucross18))] <- Ucross18$mg.u
Ucross_prep$mg.cr.u.6[which(rownames(Ucross_prep) %in% rownames(Ucross6))] <- Ucross6$mg.cr.u
Ucross_prep$mg.cr.u.12[which(rownames(Ucross_prep) %in% rownames(Ucross12))] <- Ucross12$mg.cr.u
Ucross_prep$mg.cr.u.18[which(rownames(Ucross_prep) %in% rownames(Ucross18))] <- Ucross18$mg.cr.u
Ucross_prep$ma.u.6[which(rownames(Ucross_prep) %in% rownames(Ucross6))] <- Ucross6$ma.u
Ucross_prep$ma.u.12[which(rownames(Ucross_prep) %in% rownames(Ucross12))] <- Ucross12$ma.u
Ucross_prep$ma.u.18[which(rownames(Ucross_prep) %in% rownames(Ucross18))] <- Ucross18$ma.u
Ucross_prep$ma.cr.u.6[which(rownames(Ucross_prep) %in% rownames(Ucross6))] <- Ucross6$ma.cr.u
Ucross_prep$ma.cr.u.12[which(rownames(Ucross_prep) %in% rownames(Ucross12))] <- Ucross12$ma.cr.u
Ucross_prep$ma.cr.u.18[which(rownames(Ucross_prep) %in% rownames(Ucross18))] <- Ucross18$ma.cr.u
Ucross_prep$phs.u.6[which(rownames(Ucross_prep) %in% rownames(Ucross6))] <- Ucross6$phs.u
Ucross_prep$phs.u.12[which(rownames(Ucross_prep) %in% rownames(Ucross12))] <- Ucross12$phs.u
Ucross_prep$phs.u.18[which(rownames(Ucross_prep) %in% rownames(Ucross18))] <- Ucross18$phs.u
Ucross_prep$phs.cr.u.6[which(rownames(Ucross_prep) %in% rownames(Ucross6))] <- Ucross6$phs.cr.u
Ucross_prep$phs.cr.u.12[which(rownames(Ucross_prep) %in% rownames(Ucross12))] <- Ucross12$phs.cr.u
Ucross_prep$phs.cr.u.18[which(rownames(Ucross_prep) %in% rownames(Ucross18))] <- Ucross18$phs.cr.u
