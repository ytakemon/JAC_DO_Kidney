library(tidyverse)
load("./RNAseq_data/DO188b_kidney.RData")
Upheno <- read.delim(".JAC_CS_urine_chem_v1.txt", sep = "\t")

pheno <- Upheno %>%
         mutate(
           Mouse.ID = mouse.id,
           duplicated = (duplicated(Upheno$mouse.id) | duplicated(Upheno$mouse.id, fromLast = TRUE)),
           Alb = ma.u,
           Phs = phs.u) %>%
         filter((Mouse.ID %in% annot.samples$Mouse.ID) & duplicated == FALSE) %>%
         select(Mouse.ID, Alb, Phs)
