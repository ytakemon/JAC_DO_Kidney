# Updating JAC DO Kidney ANOVA using Petr's scripts.

# ANOVA analysis using Petr's Rscirpt (takes < 30mins)
# https://github.com/simecek/TheAgingProteome/blob/master/code/anova_tests.R
cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney
# qsub -v I="input.Rdata output.csv",script=rscript_name
qsub -v I="../RNAseq_data/DO188b_kidney_noprobs.RData ../Anova_output/kidney_anova_output.csv",script=anova_tests Rsubmit_args.sh

# Calculate fdr/pval-BH for every p-value calculated from Petr's code
cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts
# qsub -v I="anova_output.csv output.csv",script=rscript_name
qsub -v I="kidney_anova_output.csv kidney_anova_fdr_output.csv",script=Kidney_ANOVA_fdr Rsubmit_args.sh

# Plot p-values and fdr from ANOVA output
cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts
qsub -v I=kidney_anova_fdr_output.csv,script=Kidney_ANOVA_pval_hist Rsubmit_args.sh

# Identify significant genes if p-value and fdr are both <= 0.05
cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts
qsub -v I="kidney_anova_fdr_output.csv kidney_anova_sig_output.csv",script=Kidney_ANOVA_significant Rsubmit_args.sh

#
