#PBS -l nodes=1:ppn=10,walltime=24:00:00

module load R/3.4.1

cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts/QTLmapping

Rscript scanBestMarker_m_Akt1_p.R $col
