#PBS -l nodes=1:ppn=1,walltime=64:00:00

module load R

cd ~/cgd/QTL_mapping/kidney_combined

Rscript scripts/qtl2blup_one.r $nparts $part
