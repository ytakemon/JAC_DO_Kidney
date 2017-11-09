setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./RNAseq_data/DO188b_kidney_noprobs.RData")

write.csv(compare, file="./QTLscan/output/eQTLintAkt1thr6_chr12.csv")
write.csv(compare, file="./QTLscan/output/eQTLintAKT1thr6_chr12.csv")
