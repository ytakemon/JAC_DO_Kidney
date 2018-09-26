# R/3.5.0 Joy in Playing
library(tidyverse)
library(readxl)
options(tibble.width = Inf)
MelkList <- "~/Desktop/Melk500sublist.xlsx"
SlopeList <- "~/Desktop/Table S4.xlsx" # our slope analysis table
MelkList <- read_excel(MelkList, sheet = "Sheet1")
SlopeList <- read_excel(SlopeList, sheet = "Change Both Significant")

# Clean up
MelkList <- MelkList %>% filter(Age == "Yes") %>% select(Symbol, Definition) %>%
  rbind(c("Rgn", "Senescence marker protein 30")) %>%
  mutate(Symbol = str_to_title(Symbol))

SlopeList[SlopeList$symbol %in% MelkList$Symbol,]
