# R/3.5.0 Joy in Playing
library(tidyverse)
library(readxl)
library(biomaRt)
options(tibble.width = Inf)

# Pick mart
mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "mmusculus_gene_ensembl",
                    verbose = TRUE)

# define file paths
RefList <- "~/Desktop/Reference_geneList.xlsx"
SlopeList <- "~/Desktop/Table S4.xlsx" # our slope analysis table

# Checkout sheet names
excel_sheets(RefList)
excel_sheets(SlopeList)

#Get mouse gene lists from refernece -------------------------------------------
mMelk <- read_excel(RefList, sheet = "Melk2005_T2_mouse") %>%
  mutate(Ref = "Melk2005") %>% distinct()
mRodwell <- read_excel(RefList, sheet = "Rodwell2004_ST3_mouse") %>%
  mutate(Ref = "Rodwell2004") %>% distinct()
hRodwell <- read_excel(RefList, sheet = "Rodwell2004_ST3_human") %>%
  mutate(Ref = "Rodwell2004") %>% distinct()
mNoris <- read_excel(RefList, sheet = "Noris2012_mouse") %>%
  mutate(Ref = "Noris2012") %>% distinct()

# compile all references
Compile <- mMelk %>% bind_rows(mRodwell) %>%
  bind_rows(mNoris) %>%
  dplyr::select("Ortholog ID", "Ortholog Gene Symbol", "Ref") %>%
  rename(EnsID = "Ortholog ID",
         Symbol = "Ortholog Gene Symbol")

# add genes for query
add_genesymbol  <- c("Rgn", "Cdkn2a", "Eif2a")
# annotate extra genes
add_genesymbol <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
              filters = c("mgi_symbol"),
              values = add_genesymbol,
              mart = mart) %>%
              mutate(Ref = "Extra") %>%
              rename(EnsID = "ensembl_gene_id",
                     Symbol = "mgi_symbol")
Compile <- bind_rows(Compile, add_genesymbol)

# Compare with slope analysis -------------------------------------------------
# read in slope dataset
Slope <- read_excel(SlopeList, sheet = "Change Both Significant") %>% arrange(symbol)
Compile <- mutate(Compile,
  INslope = EnsID %in% Slope$gene_id)

# Compare with overall ANOVA --------------------------------------------------
# read in anova dataset
mANOVA <- read.csv("~/Dropbox/TheAgingKidneyData/ANOVA/mrna.kidney_anova_table.csv")
mANOVA <- filter(mANOVA, sig.mRNA_Age.Sex == TRUE)
pANOVA <- read.csv("~/Dropbox/TheAgingKidneyData/ANOVA/protein.kidney_anova_table.csv")
pANOVA <- filter(pANOVA, sig.Prot_Age.Sex == TRUE)

# compare with anova
Compile <- mutate(Compile,
  INanova_rnasig = EnsID %in% mANOVA$id,
  INanova_protsig = EnsID %in% pANOVA$gene_id) %>%
  arrange(Symbol)

# Filter to list found in references ------------------------------------------
Gene2Include <- filter(Compile, INslope == TRUE | INanova_rnasig == TRUE | INanova_protsig == TRUE)

# in SlopeList
slopeOUT <- filter(Gene2Include, INslope == TRUE) %>% arrange(Symbol)
Slope_sub <- filter(Slope, gene_id %in% slopeOUT$EnsID)
slopeOUT <- mutate(slopeOUT,
  GroupA = Slope_sub$GroupA,
  GroupB = Slope_sub$GroupB,
  GroupC = Slope_sub$GroupC,
  GroupD = Slope_sub$GroupD,
)


write.csv(slopeOUT, file = "~/Desktop/slopeOUT.csv", row.names = FALSE)
# in RNAsig Only
RNAonlyOUT <- filter(Gene2Include, INslope == FALSE & INanova_rnasig == TRUE)
write.csv(RNAonlyOUT, file = "~/Desktop/RNAonlyOUT.csv")
# in PROTsig only
PROTonlyOUT <- filter(Gene2Include, INslope == FALSE & INanova_protsig == TRUE)
write.csv(PROTonlyOUT, file = "~/Desktop/PROTonlyOUT.csv")
