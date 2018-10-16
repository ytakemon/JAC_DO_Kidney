cd /projects/churchill-lab/projects/JAC/Takemon_DO_crosssectional_kidney/scripts/anova

# Run ANOVA and get slope of 6667 genes with both RNA and protein data
qsub -v I="../../data/Rdata/DO188b_kidney_noprobs.RData ../../results/ANOVA/kidney_anova_slope_output.csv",script=anova_tests_slope_pairs Rsubmit_args.sh

# Generate complete mRNA/Protein list (not pair subset)
# input order: Robj, mrna_output.csv, protein_output.csv
qsub -v I="DO188b_kidney_noprobs.RData m.kidney_anova_table.csv p.kidney_anova_table.csv",script=anova_tests_slope_individual Rsubmit_args.sh
