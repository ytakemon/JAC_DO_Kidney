# Check data before sending to Catherine.
Upheno <- read.csv("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Phenotype/phenotypes/JAC_DO_all_Upheno.csv")
AnnotSample <- read.csv("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Phenotype/phenotypes/JAC_DO_all_sampleinfo.csv")

identical(as.character(Upheno$Mouse.ID), as.character(AnnotSample$Mouse.ID))
