#PBS -l nodes=1:ppn=20,walltime=5:00:00

module load R/3.4.3

cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts/QTLmapping

Rscript scan1snps_diffint_m_addAkt1_m.R $col
