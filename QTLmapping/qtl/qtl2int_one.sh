#PBS -l nodes=1:ppn=1,walltime=10:00:00

module load R

cd ~/cgd/QTL_mapping/kidney_combined

Rscript scripts/qtl2int_one.r $nparts $part
