#PBS -l nodes=1:ppn=1,walltime=22:00:00

module load R

cd ~/cgd/QTL_mapping/kidney2/scripts

Rscript QTLprot.r
