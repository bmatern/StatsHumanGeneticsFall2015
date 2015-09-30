######################################################################
# R commands      Statistics for Human Genetics and Molecular Biology
# Lab 6 2015                                   University of Minnesota
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

############################
## Annotation
###########################
# source("http://www.bioconductor.org/biocLite.R)
# biocLite("BioCaseStudies)
# biocLite("Biobase")
# biocLite("annotate")
# biocLite("hgu95av2.db)

#library("BiocCaseStudies")
library("Biobase")
library("annotate")
library("hgu95av2.db")
library(ALL)
data("ALL")
bcell = grep("^B", as.character(ALL$BT))
moltyp = which(as.character(ALL$mol.biol) 
    %in% c("NEG", "BCR/ABL"))
ALL_bcrneg = ALL[, intersect(bcell, moltyp)]
ALL_bcrneg$mol.biol = factor(ALL_bcrneg$mol.biol)

data<-exprs(ALL_bcrneg)
probename<-rownames(data)
genename<-mget(probename, hgu95av2SYMBOL)
genename[1:5]
pdf("Correlation1.pdf")
plot(data[4,], data[5,], pch=16, xlab=probename[4], ylab=probename[5], main=paste("Correlation between ", probename[4], " and ", probename[5]))
dev.off()

####################
# Multiple linear regression 
####################
x<-rnorm(200)
y1<-rnorm(100, mean=1) + 2*x[1:100]
y2<-rnorm(100, mean=5) + 2*x[101:200]
ymin<-min(y1,y2) -0.1
ymax<-max(y1,y2)  + 0.1
pdf("mlr.pdf")
plot(x[1:100], y1, pch=16, ylim=c(ymin, ymax), xlab=expression(paste(X[1])), ylab="Y", main="Interaction X1X2", xlim=c(-3,2))
points(x[101:200], y2, pch=16, col=2)
abline(a=1, b=2, lty=2)
abline(a=5, b=2, col=2, lty=2)
text(0, -2.5, expression(paste("Y=1+ 2*", X[1])), col=1, bty="n", cex=1.5)
text(-1, 5, expression(paste("Y=3+2*", X[1])), col=2, bty="n", cex=1.5)
text(1.3, -3.5, expression(paste(Z, "=1")), col=1, bty="n", cex=1.5)
text(1.3, -4.3, expression(paste(Z, "=0")), col=2, bty="n", cex=1.5)
dev.off()

####################
# Interaction
####################
x<-rnorm(200)
y1<-rnorm(100, mean=1) + 2*x[1:100]
y2<-rnorm(100, mean=2) + 0.3*x[101:200]
ymin<-min(y1,y2) -0.1
ymax<-max(y1,y2)  + 0.1
pdf("interaction.pdf")
plot(x[1:100], y1, pch=16, ylim=c(ymin, ymax), xlab=expression(paste(X[1])), ylab="Y", main="Interaction X1X2", xlim=c(-4,2))
points(x[101:200], y2, pch=16, col=2)
abline(a=1, b=2, lty=2)
abline(a=2, b=0.3, col=2, lty=2)
text(-0.5, -2.5, expression(paste("Y=1+ 2*", X[1])), col=1, bty="n", cex=1.5)
text(-2.8, 2, expression(paste("Y=2+0.3*", X[1])), col=2, bty="n", cex=1.5)
text(1.3, -3.8, expression(paste(Z, "=1")), col=1, bty="n", cex=1.5)
text(1.3, -4.3, expression(paste(Z, "=0")), col=2, bty="n", cex=1.5)
dev.off()

#######################
# regression example 
########################
int<-as.numeric(ALL_bcrneg$mol.biol) * data[5,]
fit1<-lm(data[4,] ~ data[5,] + ALL_bcrneg$mol.biol + int)
fitout<-summary(fit1)

fit2<-lm(data[4,] ~ data[5,])
aa<-summary(fit2)
pdf("elm.pdf")
plot(data[5,], data[4,], pch=16, xlab=probename[5], ylab=probename[4], main=paste("linear association between ", probename[4], " and ", probename[5]), xlim=c(4,7), ylim=c(4.5,7))
abline(a=aa$coefficient[1,1], b=aa$coefficient[2,1], lty=2, col="grey")
text(5, 6.7, paste("Y=",round(aa$coefficient[1,1],2) ,"+", round(aa$coefficient[2,1],2), "X"), cex=2)
dev.off()


######################
# FAMuss
######################
setwd("/Users/yenyiho/Desktop/BioinfoClass/Notes/Lecture5")
fmsURL<-"http://people.umass.edu/foulkes/asg/data/FMS_data.txt"
fms<-read.delim(file=fmsURL, header=TRUE, sep="\t")
attach(fms)
whsnp<-which(colnames(fms)=="esr1_rs1042717")
whbmi<-which(colnames(fms)=="pre.BMI")
Geno<-esr1_rs1042717
trait<-as.numeric(pre.BMI > 25)
geno<-ifelse(Geno=="AA", 1, 0)
fit4<-glm(trait ~ geno, data=fms, family=binomial(link=logit))
summary(fit4)


############################################
# Connection between GLM and other tests
############################################
## Connection to T-Test
##
library("Biobase")
library("genefilter")
library("ALL")
data("ALL")
bcell = grep("^B", as.character(ALL$BT))
moltyp = which(as.character(ALL$mol.biol) 
    %in% c("NEG", "BCR/ABL"))
ALL_bcrneg = ALL[, intersect(bcell, moltyp)]
ALL_bcrneg$mol.biol = factor(ALL_bcrneg$mol.biol)
data<-exprs(ALL_bcrneg)
whp<-which(rownames(data)=="1636_g_at")
g1<-data[whp, ALL_bcrneg$mol.biol=="BCR/ABL"]
g2<-data[whp,ALL_bcrneg$mol.biol=="NEG"]
summt<-t.test(g1, g2,var.equal = T)

y<-c(g2,g1)
x<-c(rep(0, length(g2)), rep(1, length(g1)))
fit<-glm(y ~ x)
summary(fit)
## Connection to chisq.test
##
fmsURL<-"http://people.umass.edu/foulkes/asg/data/FMS_data.txt"
fms<-read.delim(file=fmsURL, header=TRUE, sep="\t")
attach(fms)
whsnp<-which(colnames(fms)=="esr1_rs1042717")
whbmi<-which(colnames(fms)=="pre.BMI")
Geno<-esr1_rs1042717
trait<-as.numeric(pre.BMI > 25)
geno<-ifelse(Geno=="AA", 1, 0)
tab<-table(trait, geno)
chisq.test(tab)
fit4<-glm(trait ~ geno, data=fms, family=binomial(link=logit))
summary(fit4)









