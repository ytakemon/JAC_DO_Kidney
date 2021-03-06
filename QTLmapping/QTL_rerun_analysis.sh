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

# Check and install missing packages
list.of.packages <- c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen", "ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages)
}

library(devtools)
install_github("rqtl/qtl2")
library(qtl2geno)
library(qtl2scan)
library(qtl2convert)
library(ggplot2)
# No errors? good.
q() #don't save workspace

# Submit scan all scripts ------------------------------------------------------
cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts/QTLmapping
bash scan_all_m.sh # Time ~ 4 hours
bash scan_all_p.sh # Time ~ 25 mins
bash scanint_all_m.sh # Time ~ 1 hour 20 mins
bash scanint_all_m2.sh # Time ~ 1 hour 20 mins
bash scanint_all_p.sh # Time ~ 1 hour 20 mins

# Confirm file counts
cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts/QTLmapping/QTLscan/
ls -l addscan_mrna/ | wc -l # 22244
ls -l addscan_prot/ | wc -l # 6717
ls -l intscan_mrna/Sex/ | wc -l # 22244
ls -l intscan_mrna/Age/ | wc -l # 22244
ls -l intscan_prot/Sex/ | wc -l # 6717
ls -l intscan_prot/Age/ | wc -l # 6717
# Looks good

# Compile eQTL and pQTL --------------------------------------------------------
cd /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Scripts/QTLmapping
qsub -v script=QTLmrna Rsubmit_args.sh #  took about 8 hours, ppn 10
qsub -v script=QTLprot Rsubmit_args.sh #

# Plot best eQTL and pQTL ------------------------------------------------------
qsub -v script=eQTL_IntAge_map Rsubmit_args.sh
qsub -v script=pQTL_IntAge_map Rsubmit_args.sh
