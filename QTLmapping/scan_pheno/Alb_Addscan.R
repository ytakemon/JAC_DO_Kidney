# Data: /hpcdata/gac/derived/CGD_DO_Genoprobs/MegaMUGA_hap_probs_v6.h5
# Readme: /hpcdata/gac/derived/CGD_DO_Genoprobs/MegaMUGA_hap_probs_v6.README.txt
# Samples: /hpcdata/gac/raw/JAC_DO_Phenotypes/mouse_info/DO aging study mouse info.xlsx
library(rhdf5)
h5 <- "/hpcdata/gac/derived/CGD_DO_Genoprobs/MegaMUGA_hap_probs_v6.h5"
h5ls(h5)
jac <- h5read(file = h5, name = "/JAC") #probs are naked, need to populate dimnames
dimnames(jac$probs) <- list(jac$samples, jac$founders, jac$markers)
probs <- jac$probs



# Read in the CJB data set.
grp = h5read(file = h5filename, name = "/CJB")
names(grp)
dimnames(grp$probs) = list(grp$samples, grp$founders, grp$markers)
probs = grp$probs
rm(grp)



Jan. 19, 2017

MegaMUGA_hap_probs_64Kgrid_v6.h5 contains haplotype probabilities for 5,753 DO mice run on the MegaMUGA.

These samples are genotpyed on a subset of the MegaMUGA markers, which are available
at ftp://ftp.jax.org/MUGA/MM_snps.Rdata.

The file is in HDF5 format (http://ftp.hdfgroup.org/HDF5/) and can be accessed through several
publicly available HDF5 libraries. HDF5 is a compressed file format that allows queries by

At the end of this file, I have example code to extract one data set using R.

The data are divided into eight projects. Each project is one HDF5 "group", which is similar to
a directory.

The projects are named and are listed below. Each project group contains four data sets:
probs: 3 dimensional numeric array with dimensions num_samples x num_founders x num_markers.
samples: character vector with sample IDs.
founders: character vector with founder IDs.
markers: character vector with marker IDs.

Each cell contains the founder proportion for one sample at one marker.

The columns for each sample sum to one. i.e. all(abs(colSums(probs[1,,]) - 1.0) < 1e-8) is TRUE.

The sample IDs have a prefix appended to the beginning to specify the different projects
that each set of sample is associated with.

The projects and prefixes are:

Carol Bult        CJB     285 samples
Harrison          DEH    1008 samples
JAX Aging Center  JAC    1536 samples
Karen Svenson     KLS     293 samples
Elissa Chesler    EJC     277 samples
Leah Raw Donahue  LRD     280 samples
Alan Pack         PACK    327 samples
Ellison           ellison 636 samples
Daniel Pomp       DP      605 samples
UNC               UNC     506 samples

Below is an exmaple of how to get a listing of the file contents and extract one data set.

h5filename = "/hpcdata/gac/derived/CGD_DO_Genoprobs/MegaMUGA_hap_probs_v6.h5"
library(rhdf5)
# List the contents of the file.
h5ls(h5filename)

      group     name       otype dclass              dim
0         /      CJB   H5I_GROUP
1      /CJB founders H5I_DATASET STRING                8
2      /CJB  markers H5I_DATASET STRING            68268
3      /CJB    probs H5I_DATASET  FLOAT  285 x 8 x 68268
4      /CJB  samples H5I_DATASET STRING              285
5         /      DEH   H5I_GROUP
6      /DEH founders H5I_DATASET STRING                8
7      /DEH  markers H5I_DATASET STRING            68268
8      /DEH    probs H5I_DATASET  FLOAT 1008 x 8 x 68268
9      /DEH  samples H5I_DATASET STRING             1008
10        /       DP   H5I_GROUP
11      /DP founders H5I_DATASET STRING                8
12      /DP  markers H5I_DATASET STRING            68268
13      /DP    probs H5I_DATASET  FLOAT  605 x 8 x 68268
14      /DP  samples H5I_DATASET STRING              605
15        /      EJC   H5I_GROUP
16     /EJC founders H5I_DATASET STRING                8
17     /EJC  markers H5I_DATASET STRING            68268
18     /EJC    probs H5I_DATASET  FLOAT  277 x 8 x 68268
19     /EJC  samples H5I_DATASET STRING              277
20        /      JAC   H5I_GROUP
21     /JAC founders H5I_DATASET STRING                8
22     /JAC  markers H5I_DATASET STRING            68268
23     /JAC    probs H5I_DATASET  FLOAT 1536 x 8 x 68268
24     /JAC  samples H5I_DATASET STRING             1536
25        /      KLS   H5I_GROUP
26     /KLS founders H5I_DATASET STRING                8
27     /KLS  markers H5I_DATASET STRING            68268
28     /KLS    probs H5I_DATASET  FLOAT  293 x 8 x 68268
29     /KLS  samples H5I_DATASET STRING              293
30        /      LRD   H5I_GROUP
31     /LRD founders H5I_DATASET STRING                8
32     /LRD  markers H5I_DATASET STRING            68268
33     /LRD    probs H5I_DATASET  FLOAT  280 x 8 x 68268
34     /LRD  samples H5I_DATASET STRING              280
35        /     PACK   H5I_GROUP
36    /PACK founders H5I_DATASET STRING                8
37    /PACK  markers H5I_DATASET STRING            68268
38    /PACK    probs H5I_DATASET  FLOAT  327 x 8 x 68268
39    /PACK  samples H5I_DATASET STRING              327
40        /      UNC   H5I_GROUP
41     /UNC founders H5I_DATASET STRING                8
42     /UNC  markers H5I_DATASET STRING            68268
43     /UNC    probs H5I_DATASET  FLOAT  506 x 8 x 68268
44     /UNC  samples H5I_DATASET STRING              506
45        /  ellison   H5I_GROUP
46 /ellison founders H5I_DATASET STRING                8
47 /ellison  markers H5I_DATASET STRING            68268
48 /ellison    probs H5I_DATASET  FLOAT  636 x 8 x 68268
49 /ellison  samples H5I_DATASET STRING              636

# Read in the CJB data set.
grp = h5read(file = h5filename, name = "/CJB")
names(grp)
dimnames(grp$probs) = list(grp$samples, grp$founders, grp$markers)
probs = grp$probs
rm(grp)
