# Find, clean, and combine Aging DO cross-sectional and longitudinal data
# Note: Excluding HDO data

library(rhdf5)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra)

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

# Prepare cross-sectional data
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

Ucross_prep$study <- "Cross-sectional"

# Plot data before combining ---------------------------------------------------
# Creatinine
cr6 <- ggplot() +
  geom_histogram(data = Ucross_prep, aes(Ucross_prep$cr.u.6, fill = "red", alpha = 0.2)) +
  geom_histogram(data = Ulong, aes(Ulong$X28.wk.Cr.urine..mg.dL., fill = "blue", alpha = 0.2))+
  guides(fill = "none", alpha = "none") +
  labs(title= "cr.u.6")
cr12 <- ggplot() +
  geom_histogram(data = Ucross_prep, aes(Ucross_prep$cr.u.12, fill = "red", alpha = 0.2)) +
  geom_histogram(data = Ulong, aes(Ulong$X54.wk.Cr.urine..mg.dL., fill = "blue", alpha = 0.2))+
  guides(fill = "none", alpha = "none") +
  labs(title= "cr.u.12")
cr18 <- ggplot() +
  geom_histogram(data = Ucross_prep, aes(Ucross_prep$cr.u.18, fill = "red", alpha = 0.2)) +
  geom_histogram(data = Ulong, aes(Ulong$X80.week.Cr.urine..mg.dL., fill = "blue", alpha = 0.2))+
  guides(fill = "none", alpha = "none") +
  labs(title= "cr.u.18")
plot_dir <- "/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Plot/"
pdf(file = paste0(plot_dir, "PhenoCheck_CrU.pdf"), width = 10, heigh = 5)
grid.arrange(cr6, cr12, cr18, ncol = 3)
dev.off()

# Alb
Alb6 <- ggplot() +
  geom_histogram(data = Ucross_prep, aes(Ucross_prep$ma.u.6, fill = "red", alpha = 0.2)) +
  geom_histogram(data = Ulong, aes(Ulong$X28.wk.MA.urine..mg.dL., fill = "blue", alpha = 0.2))+
  guides(fill = "none", alpha = "none") +
  labs(title= "ma.u.6")
Alb12 <- ggplot() +
  geom_histogram(data = Ucross_prep, aes(Ucross_prep$ma.u.12, fill = "red", alpha = 0.2)) +
  geom_histogram(data = Ulong, aes(Ulong$X54.wk.MA.urine..mg.dL., fill = "blue", alpha = 0.2))+
  guides(fill = "none", alpha = "none") +
  labs(title= "ma.u.12")
Alb18 <- ggplot() +
  geom_histogram(data = Ucross_prep, aes(Ucross_prep$ma.u.18, fill = "red", alpha = 0.2)) +
  geom_histogram(data = Ulong, aes(Ulong$X80.week.MA.urine..mg.dL., fill = "blue", alpha = 0.2))+
  guides(fill = "none", alpha = "none") +
  labs(title= "ma.u.18")
pdf(file = paste0(plot_dir, "PhenoCheck_AlbU.pdf"), width = 10, heigh = 5)
grid.arrange(Alb6, Alb12, Alb18, ncol = 3)
dev.off()

#Mg
Mg6 <- ggplot() +
  geom_histogram(data = Ucross_prep, aes(Ucross_prep$mg.u.6, fill = "red", alpha = 0.2)) +
  geom_histogram(data = Ulong, aes(Ulong$X28.wk.Mg.urine..mg.dL., fill = "blue", alpha = 0.2))+
  guides(fill = "none", alpha = "none") +
  labs(title= "mg.u.6")
Mg12 <- ggplot() +
  geom_histogram(data = Ucross_prep, aes(Ucross_prep$mg.u.12, fill = "red", alpha = 0.2)) +
  geom_histogram(data = Ulong, aes(Ulong$X54.wk.Mg.urine..mg.dL., fill = "blue", alpha = 0.2))+
  guides(fill = "none", alpha = "none") +
  labs(title= "mg.u.12")
Mg18 <- ggplot() +
  geom_histogram(data = Ucross_prep, aes(Ucross_prep$mg.u.18, fill = "red", alpha = 0.2)) +
  geom_histogram(data = Ulong, aes(Ulong$X80.week.Mg.urine..mg.dL., fill = "blue", alpha = 0.2))+
  guides(fill = "none", alpha = "none") +
  labs(title= "mg.u.18")
pdf(file = paste0(plot_dir, "PhenoCheck_MgU.pdf"), width = 10, heigh = 5)
grid.arrange(Mg6, Mg12, Mg18, ncol = 3)
dev.off()

#Phs
Phs6 <- ggplot() +
  geom_histogram(data = Ucross_prep, aes(Ucross_prep$phs.u.6, fill = "red", alpha = 0.2)) +
  geom_histogram(data = Ulong, aes(Ulong$X28.wk.Phs..mg.dL., fill = "blue", alpha = 0.2))+
  guides(fill = "none", alpha = "none") +
  labs(title= "phs.u.6")
Phs12 <- ggplot() +
  geom_histogram(data = Ucross_prep, aes(Ucross_prep$phs.u.12, fill = "red", alpha = 0.2)) +
  geom_histogram(data = Ulong, aes(Ulong$X54.wk.Phs..mg.dL., fill = "blue", alpha = 0.2))+
  guides(fill = "none", alpha = "none") +
  labs(title= "phs.u.12")
Phs18 <- ggplot() +
  geom_histogram(data = Ucross_prep, aes(Ucross_prep$phs.u.18, fill = "red", alpha = 0.2)) +
  geom_histogram(data = Ulong, aes(Ulong$X80.week.Phs..mg.dL., fill = "blue", alpha = 0.2))+
  guides(fill = "none", alpha = "none") +
  labs(title= "phs.u.18")
pdf(file = paste0(plot_dir, "PhenoCheck_PhsU.pdf"), width = 10, heigh = 5)
grid.arrange(Phs6, Phs12, Phs18, ncol = 3)
dev.off()

# Looks good, can combine longitutional with cross-sectional data.
# Prepare longitudinal data
Ulong_prep <- as.data.frame(Ulong[,c("Mouse.ID")])
Ulong_prep$cr.u.6 <- Ulong$X28.wk.Cr.urine..mg.dL.
Ulong_prep$cr.u.12 <- Ulong$X54.wk.Cr.urine..mg.dL.
Ulong_prep$cr.u.18 <- Ulong$X80.week.Cr.urine..mg.dL.
Ulong_prep$mg.u.6 <- Ulong$X28.wk.Mg.urine..mg.dL.
Ulong_prep$mg.u.12 <- Ulong$X54.wk.Mg.urine..mg.dL.
Ulong_prep$mg.u.18 <- Ulong$X80.week.Mg.urine..mg.dL.
Ulong_prep$mg.cr.u.6 <- Ulong$X28.wk.Mg.urine..mg.grcreat.
Ulong_prep$mg.cr.u.12 <- Ulong$X54.wk.Mg.urine..mg.grcreat.
Ulong_prep$mg.cr.u.18 <- Ulong$X80.week.Mg.urine..mg.grcreat.
Ulong_prep$ma.u.6 <- Ulong$X28.wk.MA.urine..mg.dL.
Ulong_prep$ma.u.12 <- Ulong$X54.wk.MA.urine..mg.dL.
Ulong_prep$ma.u.18 <- Ulong$X80.week.MA.urine..mg.dL.
Ulong_prep$ma.cr.u.6 <- Ulong$X28.wk.MA.urine..mg.grcreat.
Ulong_prep$ma.cr.u.12 <- Ulong$X54.wk.MA.urine..mg.grcreat.
Ulong_prep$ma.cr.u.18 <- Ulong$X80.week.MA.urine..mg.grcreat.
Ulong_prep$phs.u.6 <- Ulong$X28.wk.Phs..mg.dL.
Ulong_prep$phs.u.12 <- Ulong$X54.wk.Phs..mg.dL.
Ulong_prep$phs.u.18 <- Ulong$X80.week.Phs..mg.dL.
Ulong_prep$phs.cr.u.6 <- Ulong$X28.wk.Phs..mg.grcreat.
Ulong_prep$phs.cr.u.12 <- Ulong$X54.wk.Phs..mg.grcreat.
Ulong_prep$phs.cr.u.18 <- Ulong$X80.weelkPhs..mg.grcreat.
Ulong_prep$study <- "longitudinal"
names(Ulong_prep)[1] <- "Mouse.ID"
rownames(Ulong_prep) <- Ulong_prep$Mouse.ID

# Combine
Ucombine <- rbind(Ucross_prep, Ulong_prep)
Ucombine$duplicated <- (duplicated(Ucombine$Mouse.ID) | duplicated(Ucombine$Mouse.ID, fromLast = TRUE))
Ucombine <- Ucombine[Ucombine$duplicated == FALSE,]
rownames(Ucombine) <- Ucombine$Mouse.ID
# Make sure Ulong matches
Ulong_prep <- Ulong_prep[rownames(Ulong_prep) %in% rownames(Ucombine),]
Ulong <- Ulong[rownames(Ulong) %in% rownames(Ulong_prep),]

# Add rest of time point from longitudinal data
Ucombine$cr.u.26 <- NA
Ucombine$cr.u.33 <- NA
Ucombine$mg.u.26 <- NA
Ucombine$mg.u.33 <- NA
Ucombine$mg.cr.u.26 <- NA
Ucombine$mg.cr.u.33 <- NA
Ucombine$ma.u.26 <- NA
Ucombine$ma.u.33 <- NA
Ucombine$ma.cr.u.26 <- NA
Ucombine$ma.cr.u.33 <- NA
Ucombine$phs.u.26 <- NA
Ucombine$phs.u.33 <- NA
Ucombine$phs.cr.u.26 <- NA
Ucombine$phs.cr.u.33 <- NA
# fill
Ucombine$cr.u.26[ which(rownames(Ucombine) %in% rownames(Ulong))] <- Ulong$X106..week.Cr.urine..mg.dL.
Ucombine$cr.u.33[ which(rownames(Ucombine) %in% rownames(Ulong))] <- Ulong$X134..week.Cr.urine..mg.dL.
Ucombine$mg.u.26[ which(rownames(Ucombine) %in% rownames(Ulong))] <- Ulong$X106.week.Mg.urine..mg.dL.
Ucombine$mg.u.33[ which(rownames(Ucombine) %in% rownames(Ulong))] <- Ulong$X134.week.Mg.urine..mg.dL.
Ucombine$mg.cr.u.26[ which(rownames(Ucombine) %in% rownames(Ulong))] <- Ulong$X106.week.Mg.urine..mg.grcreat.
Ucombine$mg.cr.u.33[ which(rownames(Ucombine) %in% rownames(Ulong))] <- Ulong$X134.week.Mg.urine..mg.grcreat.
Ucombine$ma.u.26[ which(rownames(Ucombine) %in% rownames(Ulong))] <- Ulong$X106.week.MA.urine..mg.dL.
Ucombine$ma.u.33[ which(rownames(Ucombine) %in% rownames(Ulong))] <- Ulong$X134.week.MA.urine..mg.dL.
Ucombine$ma.cr.u.26[ which(rownames(Ucombine) %in% rownames(Ulong))] <- Ulong$X106.week.MA.urine..mg.grcreat.
Ucombine$ma.cr.u.33[ which(rownames(Ucombine) %in% rownames(Ulong))] <- Ulong$X134.week.MA.urine..mg.grcreat.
Ucombine$phs.u.26[ which(rownames(Ucombine) %in% rownames(Ulong))] <- Ulong$X106.week.Phs..mg.dL.
Ucombine$phs.u.33[ which(rownames(Ucombine) %in% rownames(Ulong))] <- Ulong$X134.week.Phs..mg.dL.
Ucombine$phs.cr.u.26[ which(rownames(Ucombine) %in% rownames(Ulong))] <- Ulong$X106.weelkPhs..mg.grcreat.
Ucombine$phs.cr.u.33[ which(rownames(Ucombine) %in% rownames(Ulong))] <- Ulong$X134.weelkPhs..mg.grcreat.

# Put study and duplicated at the end
Ucombine <- Ucombine[,c(1:22,25:38,23:24)]
Ucombine <- arrange(Ucombine, Mouse.ID)

# Manually check and find typo
Ucombine$Mouse.ID[5] <- "DO-1021"
Ucombine <- arrange(Ucombine, Mouse.ID)
rownames(Ucombine) <- Ucombine$Mouse.ID

# Take intersection of genoprobs and Ucombine samples
genoprobs <- sub[rownames(sub) %in% rownames(Ucombine),,]
Upheno <- Ucombine[rownames(Ucombine) %in% rownames(genoprobs),]
samples <- samples[rownames(samples) %in% rownames(Upheno),]

# Need to convert megamuga snps to positions
# Get MegaMuga snps
load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata")) #obj name: MM_snps
save(genoprobs, Upheno, samples, MM_snps, file = "/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/RNAseq_data/DO1045_kidney.Rdata")
