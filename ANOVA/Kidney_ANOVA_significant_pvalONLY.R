# Calculate fdr from p-values generated from ANOVA for both RNA and protein:
#kidney_anova_output.csv

# Usage on Cadillac
# cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts
# qsub -v I="input.csv output.csv",script=Kidney_ANOVA_significant Rsubmit_args.sh

# R/3.3.2

# Load data
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")
args <- commandArgs(trailingOnly = TRUE)  # args <- "kidney_anova_slope_output.csv"
file <- list.files(path = "./Anova_output/", pattern = paste0("^",args[[1]]), recursive = TRUE)
df <- read.csv(paste0("./Anova_output/",file[[1]]), as.is = TRUE)
output <- args[[2]] # output <- "kidney_anova_output_slope_pvalONLYsig.csv"


# List columns of interest
p_list <- c("p.mRNA_Age.Sex" ,"p.mRNA_Sex.Age","p.Prot_Age.Sex" ,"p.Prot_Sex.Age",
            "p.mRNA_Age.SexProt" ,"p.mRNA_Sex.AgeProt","p.mRNA_Prot.SexAge" ,
            "p.mRNA_Prot.Sex","p.mRNA_Prot.Age" ,"p.Prot_Age.SexmRNA",
            "p.Prot_Sex.AgemRNA" ,"p.Prot_mRNA.SexAge","p.Prot_mRNA.Sex" ,
            "p.Prot_mRNA.Age","p.mRNA_Interaction" ,"p.Prot_Interaction",
            "p.mRNA_Interaction.Prot" ,"p.Prot_Interaction.mRNA")

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
  # for each p is it significant?
  for (p in 1:length(sig_list)){
    df[r,sig_list[p]] <- df[r,p_list[p]] <= 0.05
  }
}

# Output -----------------------------------------------------------------------
# write all output
write.csv(df, paste0("./Anova_output/", output), row.names = FALSE)

# write only all non-interaction/mediation analysis values
p_list <- c("p.mRNA_Age.Sex" ,"p.mRNA_Sex.Age",
            "p.Prot_Age.Sex" ,"p.Prot_Sex.Age",
            "p.mRNA_Interaction.Prot" ,"p.Prot_Interaction.mRNA")
p_list <- which(colnames(df) %in% p_list)
sig_list <- c("sig.mRNA_Age.Sex" ,"sig.mRNA_Sex.Age",
              "sig.Prot_Age.Sex" ,"sig.Prot_Sex.Age",
              "sig.mRNA_Interaction.Prot" ,"sig.Prot_Interaction.mRNA")
sig_list <- which(colnames(df) %in% sig_list)
slope_list <- c("m.mRNA_Age.Sex" ,"m.mRNA_Sex.Age",
              "m.Prot_Age.Sex" ,"m.Prot_Sex.Age",
              "m.mRNA_Interaction.Prot" ,"m.Prot_Interaction.mRNA")
slope_list <- which(colnames(df) %in% slope_list)
df_condensed <- df[,c(1:8, p_list, slope_list, sig_list)]
write.csv(df_condensed, paste0("./Anova_output/Condensed_", output), row.names = FALSE)
