######################################################################
# R commands      Statistics for Human Genetics and Molecular Biology
# Lab 5 2015                                  University of Minnesota
######################################################################
# This file contains the R commands for the lab.
#
# Lines beginning with the symbol '#' are comments in R.  All other
# lines contain code.
#
# In R for Windows, you may wish to open this file from the menu bar
# (File:Display file); you can then copy commands into the command
# window.  (Use the mouse to highlight one or more lines; then
# right-click and select "Paste to console".)
######################################################################


###################################################
### chunk: ANOVA 
###################################################
bcell = grep("^B", as.character(ALL$BT))
moltyp = which(as.character(ALL$mol.biol) 
    %in% c("NEG", "BCR/ABL", "ALL1/AF4"))
ALL3 = ALL[, intersect(bcell, moltyp)]
ALL3$mol.biol = factor(ALL3$mol.biol)
table(ALL3$mol.biol)

all<-exprs(ALL3)
dim(all)
whs<-which(rownames(all)=="1636_g_at")
maint<-paste("Distribution of" , rownames(all)[whs], "probe by molecular subtypes")
pdf("ANOVA.pdf")
stripchart(all[whs,] ~ ALL3$mol.biol, pch=21, vertical=TRUE, method="jitter", jitter=0.2, main=maint, ylab=rownames(all)[whs], xlab="Molecular types")
abline(v=1, col="grey", lty=2)
abline(v=2, col="grey", lty=2)
abline(v=3, col="grey", lty=2)
dev.off()


##### Nonparametric Test
summary(aov(all[whs, ] ~ ALL3$mol.biol))
kruskal.test(all[whs, ], ALL3$mol.biol, na.action=na.exclude)

#########################################################
### chunk: Summarizing and presenting categorical data
#########################################################

fmsURL<-"http://people.umass.edu/foulkes/asg/data/FMS_data.txt"
fms <- read.delim(file=fmsURL, header=T, sep="\t")
attach(fms)
dim(fms)
str(fms[,1:10])
Geno<-esr1_rs1042717
trait<-as.numeric(pre.BMI > 25)
plot(Geno)
###############################
# Calculate allele frequency
###############################
GenoCount <- table(esr1_rs1042717)
GenoCount
NumbObs <- sum(!is.na(esr1_rs1042717))
GenoFreq <- as.vector(GenoCount/NumbObs)
GenoFreq
FreqA <- (2*GenoFreq[1] + GenoFreq[2])/2
FreqA
FreqG <- (GenoFreq[2] + 2*GenoFreq[3])/2
FreqG

###############################
## Caculate frequency table
################################

tab<-table(trait, Geno)
prop.table(tab, margin=1)
plot(tab)
########################
# Pearson Chi-squre Test 
# Fisher Exact Test
# Trend Test
########################

chisq.test(tab)
fisher.test(tab)



install.packages("coin")
library(coin)
attach(fms)
Geno <- esr1_rs1042717
Trait <- as.numeric(pre.BMI>25)
GenoOrd <- ordered(Geno)
test1<-independence_test(Trait~GenoOrd,teststat="quad", scores=list(GenoOrd=c(0,1,2)))
test2<-independence_test(Trait~GenoOrd,teststat="scalar", scores=list(GenoOrd=c(0,1,2)))

################
# Hardy-Weinberg Equilibrium
#################

keep<-which(!is.na(Geno) & !is.na(pre.BMI))
library(genetics)
Geno <- esr1_rs1042717[keep]
gt<-genotype(Geno, sep="")
HWE.chisq(gt)
























