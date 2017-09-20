# QTL analysis
# Regenerating qlt using Petr's codes
# Open interactive session on Cadillac

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
