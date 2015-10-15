# Example 5.1
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Example 5.1 (EM approach to haplotype frequency estimation)
attach(fms)
install.packages("haplo.stats")
library(haplo.stats)
Geno <- cbind(substr(actn3_r577x,1,1), substr(actn3_r577x,2,2),
 	substr(actn3_rs540874,1,1), substr(actn3_rs540874,2,2),
 	substr(actn3_rs1815739,1,1), substr(actn3_rs1815739,2,2),
 	substr(actn3_1671064,1,1), substr(actn3_1671064,2,2))
SNPnames <- c("actn3_r577x", "actn3_rs540874", "actn3_rs1815739",
	"actn3_1671064")
Geno.C <- Geno[Race=="Caucasian" & !is.na(Race),]
HaploEM <- haplo.em(Geno.C, locus.label=SNPnames, 
control=haplo.em.control(min.posterior=1e-4))
HaploEM
Geno.AA <- Geno[Race=="African Am" & !is.na(Race),]
HaploEM2 <- haplo.em(Geno.AA, locus.label=SNPnames, 
control=haplo.em.control(min.posterior=1e-4))
HaploEM2

# Example 5.2
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Necessary code from Example 5.1 
attach(fms)
library(haplo.stats)
Geno <- cbind(substr(actn3_r577x,1,1), substr(actn3_r577x,2,2),
 	substr(actn3_rs540874,1,1), substr(actn3_rs540874,2,2),
 	substr(actn3_rs1815739,1,1), substr(actn3_rs1815739,2,2),
 	substr(actn3_1671064,1,1), substr(actn3_1671064,2,2))
SNPnames <- c("actn3_r577x", "actn3_rs540874", "actn3_rs1815739",
	"actn3_1671064")
Geno.C <- Geno[Race=="Caucasian" & !is.na(Race),]
HaploEM <- haplo.em(Geno.C,locus.label=SNPnames, 
control=haplo.em.control(min.posterior=1e-4))
HaploEM

# Example 5.2 (Calculating posterior haplotype probabilities)
HaploEM$nreps[1:5]
HaploEM$indx.subj[1:8]
HaploEM$hap1code[1:8]
HaploEM$hap2code[1:8]
HaploEM$post[1:8]
HapProb <- HaploEM$hap.prob
HapProb
p1 <- 2*prod(HapProb[c(3,8)])
p2 <- 2*prod(HapProb[c(4,7)])
p1 / (p1+p2)
p2 / (p1+p2)

# Example 5.3
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Necessary code from Example 5.1 
attach(fms)
library(haplo.stats)
Geno <- cbind(substr(actn3_r577x,1,1), substr(actn3_r577x,2,2),
 	substr(actn3_rs540874,1,1), substr(actn3_rs540874,2,2),
 	substr(actn3_rs1815739,1,1), substr(actn3_rs1815739,2,2),
 	substr(actn3_1671064,1,1), substr(actn3_1671064,2,2))
SNPnames <- c("actn3_r577x", "actn3_rs540874", "actn3_rs1815739",
	"actn3_1671064")
Geno.C <- Geno[Race=="Caucasian" & !is.na(Race),]
HaploEM <- haplo.em(Geno.C, locus.label=SNPnames, 
control=haplo.em.control(min.posterior=1e-4))
HaploEM

Geno.AA <- Geno[Race=="African Am" & !is.na(Race),]
HaploEM2 <- haplo.em(Geno.AA, locus.label=SNPnames, 
control=haplo.em.control(min.posterior=1e-4))
HaploEM2

# New functions needed for Example 5.3

##########################################################################
# Description: This function creates a design matrix with i,j 
#		element equal to the conditional expectation 
#		of the number of copies of haplotype j for 
#		individual i based on the output from haplo.em()
# Input:	HaploEM (object resulting from haplo.em())
# Output:	XmatHap
##########################################################################
 
HapDesign <- function(HaploEM){
	Nobs <- length(unique(HaploEM$indx.subj)) # number of observations
	Nhap <- length(HaploEM$hap.prob)	# number of haplotypes
	XmatHap <- matrix(data=0,nrow=Nobs,ncol=Nhap)
	for (i in 1:Nobs){
		IDSeq <- seq(1:sum(HaploEM$nreps))[HaploEM$indx.subj==i]
		for (j in 1:length(IDSeq)){
			XmatHap[i,HaploEM$hap1code[IDSeq][j]] <- 
				XmatHap[i,HaploEM$hap1code[IDSeq][j]] + 
				HaploEM$post[IDSeq][j]
			XmatHap[i,HaploEM$hap2code[IDSeq][j]] <- 
				XmatHap[i,HaploEM$hap2code[IDSeq][j]] + 
				HaploEM$post[IDSeq][j]
			}
		}	
	return(XmatHap)
}

##########################################################################
# Description: This function creates a vector with jth element 
#		equal to the standard error of haplotype j 
#		based on the output from haplo.em()
# Input: 	HaploEM (object resulting from haplo.em())
# Output: 	HapSE
##########################################################################
 

HapFreqSE <- function(HaploEM){
	HapMat <- HapDesign(HaploEM)
	Nobs <- length(unique(HaploEM$indx.subj)) # number of observations
	Nhap <- length(HaploEM$hap.prob)	# number of haplotypes
	S.Full<-matrix(data=0, nrow=Nobs, ncol=Nhap-1)
	for(i in 1:Nobs){
		for(k in 1:(Nhap-1)){
		S.Full[i,k]<-HapMat[i,k]/HaploEM$hap.prob[k]-
			HapMat[i,Nhap]/HaploEM$hap.prob[Nhap]
		}
	}
	Score<-t(S.Full)%*%S.Full
	invScore<-solve(Score)
	HapSE<-c(sqrt(diag(invScore)), 
		sqrt(t(rep(1,Nhap-1))%*%invScore%*%rep(1,Nhap-1)))
	return(HapSE)
	}


# Example 5.3 (Testing hypotheses about haplotype frequencies within the 
EM framework)
FreqDiff <- HaploEM2$hap.prob[4] - HaploEM$hap.prob[4]
s1 <- HapFreqSE(HaploEM)[4] 
s2 <- HapFreqSE(HaploEM2)[4]
SE <- sqrt(s1^2 + s2^2)
CI <- c(FreqDiff - 1.96*SE, FreqDiff + 1.96*SE)
CI

# Example 5.4
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Necessary code from Example 5.1 
attach(fms)
library(haplo.stats)
Geno <- cbind(substr(actn3_r577x,1,1), substr(actn3_r577x,2,2),
 	substr(actn3_rs540874,1,1), substr(actn3_rs540874,2,2),
 	substr(actn3_rs1815739,1,1), substr(actn3_rs1815739,2,2),
 	substr(actn3_1671064,1,1), substr(actn3_1671064,2,2))
SNPnames <- c("actn3_r577x", "actn3_rs540874", "actn3_rs1815739",
	"actn3_1671064")
Geno.C <- Geno[Race=="Caucasian" & !is.na(Race),]
HaploEM <- haplo.em(Geno.C,locus.label=SNPnames, 
control=haplo.em.control(min.posterior=1e-4))
HaploEM

# New function needed for Example 5.4

##########################################################################
# Description: This function creates a design matrix with i,j 
#		element equal to the conditional expectation 
#		of the number of copies of haplotype j for 
#		individual i based on the output from haplo.em()
# Input:	HaploEM (object resulting from haplo.em())
# Output:	XmatHap
##########################################################################
 
HapDesign <- function(HaploEM){
	Nobs <- length(unique(HaploEM$indx.subj)) # number of observations
	Nhap <- length(HaploEM$hap.prob)	# number of haplotypes
	XmatHap <- matrix(data=0,nrow=Nobs,ncol=Nhap)
	for (i in 1:Nobs){
		IDSeq <- seq(1:sum(HaploEM$nreps))[HaploEM$indx.subj==i]
		for (j in 1:length(IDSeq)){
			XmatHap[i,HaploEM$hap1code[IDSeq][j]] <- 
				XmatHap[i,HaploEM$hap1code[IDSeq][j]] + 
				HaploEM$post[IDSeq][j]
			XmatHap[i,HaploEM$hap2code[IDSeq][j]] <- 
				XmatHap[i,HaploEM$hap2code[IDSeq][j]] + 
				HaploEM$post[IDSeq][j]
			}
		}	
	return(XmatHap)
}

# Example 5.4 (Application of haplotype trend regression (HTR)
HapMat <- HapDesign(HaploEM) 
Trait <- NDRM.CH[Race=="Caucasian" & !is.na(Race)]
mod1 <- (lm(Trait~HapMat))
mod2 <- (lm(Trait~1))
anova(mod2,mod1)

# Example 5.5
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Necessary code from Example 5.1 
attach(fms)
library(haplo.stats)
Geno <- cbind(substr(actn3_r577x,1,1), substr(actn3_r577x,2,2),
 	substr(actn3_rs540874,1,1), substr(actn3_rs540874,2,2),
 	substr(actn3_rs1815739,1,1), substr(actn3_rs1815739,2,2),
 	substr(actn3_1671064,1,1), substr(actn3_1671064,2,2))
SNPnames <- c("actn3_r577x", "actn3_rs540874", "actn3_rs1815739",
	"actn3_1671064")
Geno.C <- Geno[Race=="Caucasian" & !is.na(Race),]
HaploEM <- haplo.em(Geno.C,locus.label=SNPnames, 
control=haplo.em.control(min.posterior=1e-4))
HaploEM

# Necessary code from Example 5.4
Trait <- NDRM.CH[Race=="Caucasian" & !is.na(Race)]

# Example 5.5 (Multiple imputation for haplotype effect estimation and 
testing)
Nobs <- sum(Race=="Caucasian", na.rm=T)
Nhap <- length(HaploEM$hap.prob)
D <- 1000
Est <- rep(0,D)
SE <- rep(0,D)
for (nimput in 1:D){
	Xmat <- matrix(data=0,nrow=Nobs,ncol=Nhap)
	for (i in 1:Nobs){
		IDSeq <- seq(1:sum(HaploEM$nreps))[HaploEM$indx.subj==i]
		if (length(IDSeq)>1){Samp <- sample(IDSeq,size=1,
			prob=HaploEM$post[IDSeq])}
		if (length(IDSeq)==1){Samp <- IDSeq}	
		Xmat[i,HaploEM$hap1code[Samp]] <-1
		Xmat[i,HaploEM$hap2code[Samp]] <-1
		}	
	h8 <- Xmat[,8]>=1
	Est[nimput] <- summary(lm(Trait~h8))$coefficients[2,1]
	SE[nimput] <- summary(lm(Trait~h8))$coefficients[2,2]
}
MeanEst <- mean(Est)
Wd <- mean(SE^2)
Bd <- (1/(D-1))*sum((Est-MeanEst)^2)
Td <- Wd + ((D+1)/D)*Bd
nu <- (D-1)*(1 + (1/(D+1))*(Wd/Bd))^2
1-pt(MeanEst/sqrt(Td),df=nu)

# Example 5.6
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Necessary code from Example 5.1
attach(fms)
Geno <- cbind(substr(actn3_r577x,1,1), substr(actn3_r577x,2,2),
 	substr(actn3_rs540874,1,1), substr(actn3_rs540874,2,2),
 	substr(actn3_rs1815739,1,1), substr(actn3_rs1815739,2,2),
 	substr(actn3_1671064,1,1), substr(actn3_1671064,2,2))
Geno.C <- Geno[Race=="Caucasian" & !is.na(Race),]

# Example 5.6 (EM for estimation and testing of haplotype-trait 
association)
library(haplo.stats)
Geno.C <- setupGeno(Geno.C)
Trait <- NDRM.CH[Race=="Caucasian" & !is.na(Race)]
Dat <- data.frame(Geno.C=Geno.C, Trait=Trait)
haplo.glm(Trait~Geno.C,data=Dat, 
allele.lev=attributes(Geno.C)$unique.alleles)
haplo.glm(Trait~Geno.C,data=Dat, 
allele.lev=attributes(Geno.C)$unique.alleles, 
	control=haplo.glm.control(haplo.base=9))
haplo.glm(Trait~Geno.C,data=Dat, 
allele.lev=attributes(Geno.C)$unique.alleles,
 	control=haplo.glm.control(haplo.effect="dominant"))


