#PBS -l nodes=1:ppn=10,walltime=12:00:00

module load R/3.4.4
cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts/QTLmapping/pbatch_qc

Rscript scan_one_p.R $col
