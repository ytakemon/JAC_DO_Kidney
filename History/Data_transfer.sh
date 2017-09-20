# Copy data to my directory on Cadillac ----------------------------------------

#RNAseq data
cd /projects/cgd/QTL_mapping/kidney2/
cp DO188b_kidney* /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/RNAseq_data

# Annotation of data contained in RData:
# DO188b_kidney.RData:
## G - 188x188 matrix ?????? kinship?
## Glist - G by chromosomes
## N - list containing number of attributes
## annot.mrna - df mrna annotation (gene id, chr, gene sybol, etc...)
## annot.protein - df protein annoation
## annot.samples - df sample annotation
## covar - 1/0 covariate list
## expr.mrna - mRNA counts are they? how have they been normalized?
## expr.protein - same questions as above
## genoprobs - Genotype probability of all 188 animals
## raw.mrna - Raw mRNA reads (raw... counts?)
## raw.protein - Raw protein
## snps - 64K Marker id, chrom, posisiton

# phenotype data
cd /hpcdata/gac/projects/JAC_DO_Phenotypes
cp -r * /projects/korstanje-lab/ytakemon/JAC_DO_Kidney/Phenotype/

# Copy Petr's codes from cadillac to desktop GitHub ----------------------------

# Copy QTL mapping codes
cd /data/cgd/QTL_mapping/kidney2
scp -r scripts/ ytakemon@mdg-gandalf:~/GitHub/JAC_DO_Kidney/Petrs_codes/QTLmapping
