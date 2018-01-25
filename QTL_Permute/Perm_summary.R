setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney")

# mRNA
add1 <- readRDS("./QTLperm/eQTL_Napsa_1000perm.rds")
add2 <- readRDS("./QTLperm/eQTL_Sf3a1_1000perm.rds")
summary(add1)
summary(add2)

int1 <- readRDS("./QTLperm/eQTL_Ageint_Napsa_1000perm.rds")
int2 <- readRDS("./QTLperm/eQTL_Ageint_Sf3a1_1000perm.rds")
summary(int1)
summary(int2)

lodDiff <- data.frame(perm = c("Napsa", "Sf3a1"),
                      add = c(add1[1], add2[1]),
                      int = c(int1[1], int2[1]),
                      diff = c((int1[1] - add1[1]), (int2[1] - add2[1])))
#protein
add1 <- readRDS("./QTLperm/pQTL_Inmt_1000perm.rds")
add2 <- readRDS("./QTLperm/pQTL_Paf1_1000perm.rds")
summary(add1)
summary(add2)

int1 <- readRDS("./QTLperm/pQTL_Ageint_Inmt_1000perm.rds")
int2 <- readRDS("./QTLperm/pQTL_Ageint_Paf1_1000perm.rds")
summary(int1)
summary(int2)


lodDiff2 <- data.frame(perm = c("Inmt", "Paf1"),
                      add = c(add1[1], add2[1]),
                      int = c(int1[1], int2[1]),
                      diff = c((int1[1] - add1[1]), (int2[1] - add2[1])))

# summary is esentially a quantile(lodscores,0.95)

> lodDiff
   perm      add       int     diff
1 Napsa 7.769087 10.577989 2.808902
2 Sf3a1 7.393432  9.696623 2.303191

> lodDiff2
  perm      add      int       diff
1 Inmt 6.330818 9.006886  2.6760673
2 Paf1 7.594473 6.937915 -0.6565582
