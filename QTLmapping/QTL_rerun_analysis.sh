# QTL analysis
# Regenerating qlt using Petr's codes
# Open interactive session on Cadillac

# Create batch qsub scripts ----------------------------------------------------
# load python to run .py scripts
module load python/2.7.3

cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts/QTLmapping
# run .py scripts to generate batch bash scripts
python scan_all_m.py
python scan_all_p.py
python scanint_all_m.py
python scanint_all_m2.py
python scanint_all_p.py

# This will create the following bash scripts:
# scan_all_m.sh
# scan_all_p.sh
# scanint_all_m.sh
# scanint_all_m2.sh
# scanint_all_p.sh

# For files scan*one.R, change path to fit my directories
# For files scan*one.sh, change path of rscript and rename .r to .R

# Preparing R for usage --------------------------------------------------------
# Install packages in R/3.4.1 to use qtl2 packages
# Begin interactive session on Cadillac
module load R/3.4.1
R # R console

# Install packages
install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen"))
library(devtools)
install_github("rqtl/qtl2")
# load libraries to test it works
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
# No error? good.



bash scan_all_m.sh
