######################################################################
# R commands      Statistics for Human Genetics and Molecular Biology
# Lab 4 2015                                  University of Minnesota
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
source("http://www.bioconductor.org/biocLite.R")
biocLite("Biobase")
biocLite("genefilter")
biocLite("ALL")

library("Biobase")
library("genefilter")
library("ALL")
data("ALL")

#############################################
## Identifying tumor molecular subgroups
#############################################
bcell = grep("^B", as.character(ALL$BT))
moltyp = which(as.character(ALL$mol.biol) 
    %in% c("NEG", "BCR/ABL"))
ALL_bcrneg = ALL[, intersect(bcell, moltyp)]
ALL_bcrneg$mol.biol = factor(ALL_bcrneg$mol.biol)
data<-exprs(ALL_bcrneg)
whp<-which(rownames(data)=="1636_g_at")

###################################################
### chunk: Two Sample T-Test
###################################################
g1<-data[whp, ALL_bcrneg$mol.biol=="BCR/ABL"]
g2<-data[whp,ALL_bcrneg$mol.biol=="NEG"]
summt<-t.test(g1, g2)
tobs<-summt$statistic


#######################
# Wilcoxon Test
#######################
wilcox.test(g1,g2)


###################################################
### Exercise: Generate permuted null distribution 
###################################################



 


###################################################
### chunk: stripchart
###################################################

whp<-which(rownames(data)=="1636_g_at")
maint<-paste("Distribution of" , rownames(data)[whp], "probe by cancer molecular subtypes")
pdf("GeneDot.pdf")
stripchart(data[whp,] ~ ALL_bcrneg$mol.biol, pch=21, method="jitter", jitter=0.2, vertical=TRUE, ylab=rownames(data)[whp], main=maint, xlab="Molecular types")
abline(v=1, col="grey", lty=2)
abline(v=2, col="grey", lty=2)
lines(x=c(0.9, 1.1), rep(mean(data[whp,ALL_bcrneg$mol.biol=="BCR/ABL"]),2), col=4)
lines(x=c(1.9, 2.1), rep(mean(data[whp,ALL_bcrneg$mol.biol=="NEG"]),2), col=4)
dev.off()


###################################################
### chunk: Men & Women example
###################################################

setwd("/Users/yenyiho/Desktop/BioinfoClass/Notes/Lecture4")
f<-rnorm(15, mean=65, sd=2)
m<-rnorm(15, mean=70, sd=4)

data<-as.data.frame(cbind(c(f,m), rep(c("Female", "Male"), each=15)))
data[,1]<-as.numeric(data[,1])
pdf("sampleMen.pdf")
stripchart(data[,1] ~ data[,2], vertical=TRUE, pch=16, col=c(2,4), method="jitter", jitter=0.1, ylab="Height", main="Sample of 15 women and 15 men", cex.main=2)
abline(v=1, lty=2, col=2)
abline(v=2, lty=2, col=4)
lines(x=c(0.9, 1.1), rep(mean(data[1:15,1]),2), col=2)
lines(x=c(1.9, 2.1), rep(mean(data[16:30,1]),2), col=4)
dev.off()




