#!/bin/bash -l
#PBS -l nodes=1:ppn=3,walltime=03:00:00

# Usage:
# qsub -v I=${from list},script=${R script name} Rsubmit_args.sh

cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts
module load R/3.4.1

# I is the input variable
R --no-save --args ${I} < ${script}.R > ${script}.Rout
