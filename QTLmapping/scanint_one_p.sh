#PBS -l nodes=1:ppn=10,walltime=12:00:00

module load Rscript scanint_one_m.R $col

cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts/QTLmapping

Rscript scanint_one_p.R $col
