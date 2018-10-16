library(tidyverse)
# Extract raw RNA counts and normalized counts for GEO deposit
load("./DO188b_kidney.RData")

GEO_mrna_raw <- rownames_to_column(as.data.frame(raw.mrna)) %>%
  rename(SampleID = rowname)

write.csv(GEO_mrna_raw, "./GEO_submission/GEO_mrna_raw.csv",
  quote = FALSE, row.names = FALSE)

GEO_mrna_rank <- rownames_to_column(as.data.frame(expr.mrna)) %>%
  rename(SampleID = rowname)

write.csv(GEO_mrna_rank, "./GEO_submission/GEO_mrna_ranked.csv",
  quote = FALSE, row.names = FALSE)

# Extract lod peak data
rm(list = ls())
load("./DO188b_kidney_20180926_YT.Rdata")
