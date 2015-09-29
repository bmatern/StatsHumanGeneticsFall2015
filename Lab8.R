# Example 3.1
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Example 3.1 (Measuring LD using D-prime):
library(genetics)
attach(fms)
actn3_r577x[1:10]
actn3_rs540874[1:10]
Actn3Snp1 <- genotype(actn3_r577x,sep="")
Actn3Snp2 <- genotype(actn3_rs540874,sep="")
Actn3Snp1[1:10]
LD(Actn3Snp1,Actn3Snp2)$"D'"
Esr1Snp1 <- genotype(esr1_rs1801132,sep="")
LD(Actn3Snp1,Esr1Snp1)$"D'"

#Example 3.2
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Necessary code from Example 3.1:
attach(fms)
library(genetics)
Actn3Snp1 <- genotype(actn3_r577x,sep="")
Actn3Snp2 <- genotype(actn3_rs540874,sep="")

# Example 3.2 (Measuring pairwise LD for a group of SNPs):
Actn3Snp3 <- genotype(actn3_rs1815739,sep="")
Actn3Snp4 <- genotype(actn3_1671064,sep="")
Actn3AllSnps <- data.frame(Actn3Snp1,Actn3Snp2,Actn3Snp3,Actn3Snp4)
LD(Actn3AllSnps)$"D'"
install.packages("LDheatmap")
library(LDheatmap)
LDheatmap(Actn3AllSnps, LDmeasure="D'")

#Example 3.3
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Example 3.3 (Measuring LD based on r^2 and the \chi^2-statistic):
attach(fms)
library(genetics)
Actn3Snp1 <- genotype(actn3_r577x,sep="")
Actn3Snp2 <- genotype(actn3_rs540874,sep="")
LD(Actn3Snp1,Actn3Snp2)$"R^2"
LD(Actn3Snp1,Actn3Snp2)

# Example 3.4
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Necessary code from Example 3.2:
attach(fms)
library(genetics)
Actn3Snp1 <- genotype(actn3_r577x,sep="")
Actn3Snp2 <- genotype(actn3_rs540874,sep="")
Actn3Snp3 <- genotype(actn3_rs1815739,sep="")
Actn3Snp4 <- genotype(actn3_1671064,sep="")
Actn3AllSnps <- data.frame(Actn3Snp1,Actn3Snp2,Actn3Snp3,Actn3Snp4)

# Example 3.4 (Determining average LD across multiple SNPs):
LDMat <- LD(Actn3AllSnps)$"D'"
mean(LDMat,na.rm=T)

# Example 3.5
# Example 3.5 (Population substructure and LD):
ObsCount <- matrix(c(136,64,64,136),2)
ObsCount
ExpCount <- chisq.test(ObsCount)$expected
ExpCount

# Example 3.6
# Reading in HGDP data:
hgdp <- 
read.delim("http://people.umass.edu/foulkes/asg/data/HGDP_AKT1.txt", 
header=T, sep="\t")

# Example 3.6 (Testing for HWE using Pearsons \chi^2-test):
attach(hgdp)
Akt1Snp1 <- AKT1.C0756A
ObsCount <- table(Akt1Snp1)
Nobs <- sum(ObsCount)
ObsCount
FreqC <- (2 * ObsCount[3] + ObsCount[2])/(2*Nobs)
ExpCount <- c(Nobs*(1-FreqC)^2, 2*Nobs*FreqC*(1-FreqC),Nobs*FreqC^2)
ExpCount
ChiSqStat <- sum((ObsCount - ExpCount)^2/ExpCount)
ChiSqStat
qchisq(1-0.05,df=1)

library(genetics)
Akt1Snp1 <- genotype(AKT1.C0756A, sep="")
HWE.chisq(Akt1Snp1)

# Example 3.7
# Reading in HGDP data:
hgdp <- 
read.delim("http://people.umass.edu/foulkes/asg/data/HGDP_AKT1.txt", 
header=T, sep="\t")

# Example 3.7 (Testing for HWE using Fishers exact test):
attach(hgdp)
Akt1Snp1Maya <- AKT1.C0756A[Population=="Maya"]
ObsCount <- table(Akt1Snp1Maya)
ObsCount
Nobs <- sum(ObsCount)
FreqC <- (2 * ObsCount[3] + ObsCount[2])/(2*Nobs)
ExpCount <- c(Nobs*(1-FreqC)^2, 2*Nobs*FreqC*(1-FreqC),Nobs*FreqC^2)
ExpCount
n11 <- ObsCount[3]
n12 <- ObsCount[2]
n22 <- ObsCount[1]
n1 <- 2*n11+n12
Num <- 2^n12 * factorial(Nobs)/prod(factorial(ObsCount))
Denom <- factorial(2*Nobs) / (factorial(n1)*factorial(2*Nobs-n1))
FisherP1 <- Num/Denom
FisherP1

library(genetics)
Akt1Snp1Maya <- genotype(AKT1.C0756A[Population=="Maya"], sep="")
HWE.exact(Akt1Snp1Maya)

# Example 3.8
# Reading in HGDP data:
hgdp <- 
read.delim("http://people.umass.edu/foulkes/asg/data/HGDP_AKT1.txt", 
header=T, sep="\t")

# Example 3.8 (HWE and geographic origin):
attach(hgdp)
table(Geographic.area)
library(genetics)
Akt1Snp1 <- genotype(AKT1.C0756A, sep="")
HWEGeoArea <- tapply(Akt1Snp1,INDEX=Geographic.area,HWE.chisq)  
HWEGeoArea$"Central Africa"
HWEGeoArea$"South America"

#Example 3.9
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Example 3.9 (Generating a similarity matrix):
attach(fms)
NamesAkt1Snps <- names(fms)[substr(names(fms),1,4)=="akt1"]
NamesAkt1Snps
FMSgeno <- fms[,is.element(names(fms),NamesAkt1Snps)]
FMSgenoNum <- data.matrix(FMSgeno)
FMSgenoNum[is.na(FMSgenoNum)] <- 4
DistFmsGeno <- as.matrix(dist(FMSgenoNum))
DistFmsGeno[1:5,1:5]

#Example 3.10
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Necessary code from Example 3.9:
attach(fms)
NamesAkt1Snps <- names(fms)[substr(names(fms),1,4)=="akt1"]
FMSgeno <- fms[,is.element(names(fms),NamesAkt1Snps)]
FMSgenoNum <- data.matrix(FMSgeno)
FMSgenoNum[is.na(FMSgenoNum)] <- 4
DistFmsGeno <- as.matrix(dist(FMSgenoNum))

# Exammple 3.10 (Multidimensional scaling (MDS) for identifying population 
substructure):
plot(cmdscale(DistFmsGeno),xlab="C1",ylab="C2")
abline(v=0,lty=2)	
abline(h=4,lty=2)

# Example 3.11
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Necessary code from Example 3.9:
attach(fms)
NamesAkt1Snps <- names(fms)[substr(names(fms),1,4)=="akt1"]
FMSgeno <- fms[,is.element(names(fms),NamesAkt1Snps)]
FMSgenoNum <- data.matrix(FMSgeno)
FMSgenoNum[is.na(FMSgenoNum)] <- 4

# Exammple 3.11 (Principal components analysis (PCA) for identifying 
population substructure):
PCFMS <- prcomp(FMSgenoNum)
plot(PCFMS$"x"[,1],PCFMS$"x"[,2],xlab="PC1",ylab="PC2")

