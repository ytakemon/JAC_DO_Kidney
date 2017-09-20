#PBS -l nodes=1:ppn=10,walltime=12:00:00

module load R

cd ~/cgd/QTL_mapping/kidney2

Rscript scripts/scanint_one_m.r $col
