cd /projects/churchill-lab/projects/JAC/Takemon_DO_crosssectional_kidney/scripts/anova

qsub -v I="../../data/Rdata/DO188b_kidney_noprobs.RData ../../results/ANOVA/kidney_anova_slope_output.csv",script=anova_tests_slope_pairs Rsubmit_args.sh
