#PBS -l nodes=1:ppn=1,walltime=12:00:00

module load R

cd ~/cgd/QTL_mapping/kidney2

Rscript scripts/scan_one_p.r $col
