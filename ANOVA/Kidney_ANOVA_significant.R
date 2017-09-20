# Calculate fdr from p-values generated from ANOVA for both RNA and protein:
#kidney_anova_output.csv

# Usage on Cadillac
# cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts
# qsub -v I="input.csv output.csv",script=Kidney_ANOVA_significant Rsubmit_args.sh

# R/3.3.2

# Load data
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
args <- commandArgs(trailingOnly = TRUE)  # args <- kidney_anova_fdr_output.csv""
file <- list.files(pattern = paste0("^",args[[1]]), recursive = TRUE)
df <- read.csv(file[[1]], as.is = TRUE)
output <- args[[2]]


# List columns of interest
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

sig_list <- c("sig.mRNA_Age.Sex" ,"sig.mRNA_Sex.Age","sig.Prot_Age.Sex" ,"sig.Prot_Sex.Age",
            "sig.mRNA_Age.SexProt" ,"sig.mRNA_Sex.AgeProt","sig.mRNA_Prot.SexAge" ,
            "sig.mRNA_Prot.Sex","sig.mRNA_Prot.Age" ,"sig.Prot_Age.SexmRNA",
            "sig.Prot_Sex.AgemRNA" ,"sig.Prot_mRNA.SexAge","sig.Prot_mRNA.Sex" ,
            "sig.Prot_mRNA.Age","sig.mRNA_Interaction" ,"sig.Prot_Interaction",
            "sig.mRNA_Interaction.Prot" ,"sig.Prot_Interaction.mRNA")

# Add new columns / rest columns
df[, sig_list] <- NA
# Identify significant data for ever row
for (r in 1:dim(df)[1]){
  # for each pair of p and fdr values
  for (p in 1:length(sig_list)){
    df[r,sig_list[p]] <- df[r,p_list[p]] && df[r,fdr_list[p]] <= 0.05
  }
}

write.csv(df, paste0("./Anova_output/", output), row.names = FALSE)
