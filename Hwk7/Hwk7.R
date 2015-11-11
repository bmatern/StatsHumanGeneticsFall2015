# 1) 

source("https://bioconductor.org/biocLite.R")
biocLite()
#biocLite(c("GenomicFeatures", "AnnotationDbi"))
#biocLite(c("affy", "oligo", "limma", "siggenes", "pd.hg.u133a","pd.hg.u133.plus.2."))
biocLite(c("hgu133plus2.db","hgu133a.db"))
getwd()
setwd("/Users/bmatern/school/Fall2015/Stats4HumGenomics/StatsHumanGeneticsFall2015/Hwk7/GSE18088_RAW")

library(oligo)
library(siggenes)
library(limma)
library(pd.hg.u133.plus.2)
library(hgu133plus2.db)
library(hgu133a.db)

CELfiles=dir(pattern=".gz$")
rawdata=read.celfiles(filenames= CELfiles)

pdf("/Users/bmatern/school/Fall2015/Stats4HumGenomics/StatsHumanGeneticsFall2015/Hwk7/MAPlotBeforeNorm.pdf")
MAplot(rawdata)
dev.off()

dim(rawdata)
#pdf("/Users/bmatern/school/Fall2015/Stats4HumGenomics/StatsHumanGeneticsFall2015/Hwk7/Images.pdf")
#hist(rawdata)
#for (i in 1:ncol(rawdata)){
#  image(rawdata[,i])
#}

pdf("/Users/bmatern/school/Fall2015/Stats4HumGenomics/StatsHumanGeneticsFall2015/Hwk7/BoxplotBeforeNorm.pdf")
boxplot(rawdata)
dev.off()


normdata = rma(rawdata)
#getwd()
#getwd()

pdf("/Users/bmatern/school/Fall2015/Stats4HumGenomics/StatsHumanGeneticsFall2015/Hwk7/MAafterNorm.pdf")
MAplot(normdata)
dev.off()

pdf("/Users/bmatern/school/Fall2015/Stats4HumGenomics/StatsHumanGeneticsFall2015/Hwk7/BoxplotAfterNorm.pdf")
boxplot(normdata)
dev.off()
