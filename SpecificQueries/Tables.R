load("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/RNAseq_data/DO188b_kidney.RData")

f <- annot.samples[annot.samples$Sex == "F",]
f6 <- f[f$Age == 6,]
f12 <- f[f$Age == 12,]
f18 <- f[f$Age == 18,]
table(f6$Generation)
table(f12$Generation)
table(f18$Generation)
m <- annot.samples[annot.samples$Sex == "M",]
m6 <- m[m$Age == 6,]
m12 <- m[m$Age == 12,]
m18 <- m[m$Age == 18,]
table(m6$Generation)
table(m12$Generation)
table(m18$Generation)
