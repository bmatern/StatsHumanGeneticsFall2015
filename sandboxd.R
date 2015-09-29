# Example 1.1
# Reading in FAMuSS data
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
                  header=T, sep="\t")

names(fms)
names(fms)[1:20]

table(fms$Gender)
attach(fms)
table(Gender)

# Example 1.1 (Identifying the minor allele and its frequency):
GenoCount <- table(actn3_rs540874)
GenoCount
NumbObs <- sum(!is.na(actn3_rs540874))
GenoFreq <- as.vector(GenoCount/NumbObs)
GenoFreq
FreqA <- (2*GenoFreq[1] + GenoFreq[2])/2
FreqA
FreqG <- (GenoFreq[2] + 2*GenoFreq[3])/2
FreqG

# Alternatively, using genetics package:
install.packages("genetics")
library(genetics) 
Geno <- genotype(actn3_rs540874,sep="")
summary(Geno)

hgdpURL <- "http://people.umass.edu/foulkes/asg/data/HGDP_AKT1.txt"
hgdp <- read.delim(file=hgdpURL)
head(hgdp)
dim(hgdp)
length(table(hgdp[,4]))

vircoURL <- "http://people.umass.edu/foulkes/asg/data/Virco_data.csv"
virco <- read.csv(file=vircoURL)
dim(virco)
virco[1:5,c(1,6,11,32,85,93,104,112,122)]