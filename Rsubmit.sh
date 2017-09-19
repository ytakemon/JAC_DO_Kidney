#!/bin/bash -l
#PBS -l nodes=1:ppn=3,walltime=03:00:00

# Usage:
# qsub -v script=${R script name without ".R"} Rsubmit.R

cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts
module load R/3.3.2

# $script  is the input R script name
R --no-save < ${script}.R > ${script}.Rout
