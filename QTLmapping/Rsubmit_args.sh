#!/bin/bash -l
#PBS -l nodes=1:ppn=10,walltime=24:00:00

# Usage:
# qsub -v I=${from list},script=${R script name} Rsubmit_args.sh

cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts/QTLmapping
module load R/3.4.1

# I is the input variable
R --no-save --args ${I} < ${script}.R > ${script}.Rout
