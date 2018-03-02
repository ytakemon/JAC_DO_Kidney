# Updating JAC DO Kidney ANOVA using Petr's scripts.

# ANOVA: Changes in mRNA/Protein with Age / Sex as linear ----------------------
# ANOVA analysis using Petr's Rscirpt (takes < 30mins)
# https://github.com/simecek/TheAgingProteome/blob/master/code/anova_tests.R
cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts/ANOVA
# qsub -v I="input.Rdata output.csv",script=rscript_name
qsub -v I="../../RNAseq_data/DO188b_kidney_noprobs.RData ../../Anova_output/kidney_anova_output.csv",script=ANOVA/anova_tests_pairs Rsubmit_args.sh
# Calculate fdr/pval-BH for every p-value calculated from Petr's anova code
# qsub -v I="anova_output.csv output.csv",script=rscript_name
qsub -v I="kidney_anova_output.csv kidney_anova_fdr_output.csv",script=Kidney_ANOVA_fdr Rsubmit_args.sh
# Plot p-values and fdr from ANOVA output
qsub -v I=kidney_anova_fdr_output.csv,script=Kidney_ANOVA_pval_fdr_hist Rsubmit_args.sh
# Identify significant genes if p-value and fdr are both <= 0.05
qsub -v I="kidney_anova_fdr_output.csv kidney_anova_sig_output.csv",script=Kidney_ANOVA_significant_pval_FDR Rsubmit_args.sh

# Modified ANOVA analysis to inlude slope of regression line.
qsub -v I="../../RNAseq_data/DO188b_kidney_noprobs.RData ../../Anova_output/kidney_anova_slope_output.csv",script=anova_tests_slope_pairs Rsubmit_args.sh
# Plot slopes
qsub -v I=kidney_anova_slope_output.csv,script=ANOVA_slope_plot_pairs Rsubmit_args.sh

# Generate complete mRNA/Protein list (not pair subset)
# input order: Robj, mrna_output.csv, protein_output.csv
qsub -v I="DO188b_kidney_noprobs.RData m.kidney_anova_table.csv p.kidney_anova_table.csv",script=anova_tests_slope_pairs Rsubmit_args.sh
# Plot slopes
qsub -v I="mrna.kidney_anova_table.csv protein.kidney_anova_table.csv",script=ANOVA_slope_plot_individual Rsubmit_args.sh





# Sharing gene lists with Gary & Ron -------------------------------------------
cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/ANOVA
# qsub -v I="input.Rdata output.csv",script=rscript_name
qsub -v I="../../RNAseq_data/DO188b_kidney_noprobs.RData ../../Anova_output/kidney_anova_output.csv",script=anova_tests_pairs Rsubmit_args.sh
# Identify significant genes if p-value and fdr are both <= 0.05
qsub -v I="kidney_anova_fdr_output.csv kidney_anova_sig_output.csv",script=Kidney_ANOVA_significant_pvalONLY Rsubmit_args.sh
# Plot slopes, includes gene list output
qsub -v I=kidney_anova_slope_output.csv,script=ANOVA_slope_plot Rsubmit_args.sh














# FACTOR ANALYSIS DISCONTINUED
# ANOVA: Changes in mRNA/Protein with Age as factor ----------------------------
cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/ANOVA
# ANOVA mRNA/Protein with Age as factor. Generates p-values and slopes
qsub -v I="../../RNAseq_data/DO188b_kidney_noprobs.RData ../../Anova_output/kidney_anova_output.csv",script=anova_tests_asfactor Rsubmit_args.sh
# Plot pval histgrams
qsub -v I=kidney_anova_factor_table.csv,script=Kidney_ANOVA_pval_hist_asfactor Rsubmit_args.sh
