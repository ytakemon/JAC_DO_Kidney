# Open interactive session on Cadillac
cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney

# Load modules
# module load R/3.3.2

# ANOVA analysis using Petr's Rscirpt (take < 30mins)
# https://github.com/simecek/TheAgingProteome/blob/master/code/anova_tests.R
Rscript ./Scripts/anova_tests.R ./RNAseq_data/DO188b_kidney_noprobs.RData ./Anova_output/kidney_anova_output.csv

# Plot ANOVA output
cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts
qsub -v I=kidney_anova_output.csv,script=Kidney_ANOVA_pval_hist Rsubmit_args.sh
