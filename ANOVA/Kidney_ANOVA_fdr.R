# Calculate fdr from p-values generated from ANOVA for both RNA and protein:
#kidney_anova_output.csv

# Usage on Cadillac
# cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts
# qsub -v I=kidney_anova_output.csv,script=Kidney_ANOVA_fdr Rsubmit_args.R

# R/3.3.2

# Load data
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
args <- commandArgs(trailingOnly = TRUE)
anova <- list.files(pattern = paste0("^",args[[1]]), recursive = TRUE)
anova <- read.csv(anova[[1]])
output <- args[[2]]


p_list <- c("p.mRNA_Age.Sex" ,"p.mRNA_Sex.Age","p.Prot_Age.Sex" ,"p.Prot_Sex.Age",
            "p.mRNA_Age.SexProt" ,"p.mRNA_Sex.AgeProt","p.mRNA_Prot.SexAge" ,
            "p.mRNA_Prot.Sex","p.mRNA_Prot.Age" ,"p.Prot_Age.SexmRNA",
            "p.Prot_Sex.AgemRNA" ,"p.Prot_mRNA.SexAge","p.Prot_mRNA.Sex" ,
            "p.Prot_mRNA.Age","p.mRNA_Interaction" ,"p.Prot_Interaction",
            "p.mRNA_Interaction.Prot" ,"p.Prot_Interaction.mRNA")

fdr_list <- c("fdr.mRNA_Age.Sex" ,"fdr.mRNA_Sex.Age","fdr.Prot_Age.Sex" ,"fdr.Prot_Sex.Age",
            "fdr.mRNA_Age.SexProt" ,"fdr.mRNA_Sex.AgeProt","fdr.mRNA_Prot.SexAge" ,
            "fdr.mRNA_Prot.Sex","fdr.mRNA_Prot.Age" ,"fdr.Prot_Age.SexmRNA",
            "fdr.Prot_Sex.AgemRNA" ,"fdr.Prot_mRNA.SexAge","fdr.Prot_mRNA.Sex" ,
            "fdr.Prot_mRNA.Age","fdr.mRNA_Interaction" ,"fdr.Prot_Interaction",
            "fdr.mRNA_Interaction.Prot" ,"fdr.Prot_Interaction.mRNA")

for (i in 1:length(p_list)){
  anova[, fdr_list[i]] <- p.adjust( anova[, p_list[i]],
                                     method = "BH",
                                     n = length(anova[, p_list[i]]))
}

write.csv(anova, paste0("./Anova_output/", output), row.names = FALSE)
