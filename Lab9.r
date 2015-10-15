# Example 4.1
# Reading in Virco data:
virco <- 
read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", 
header=T, sep=",")

# Example 4.1 (Bonferroni adjustment):
attach(virco)
PrMut <- virco[,23:121]!="-" & virco[,23:121]!="."
NObs <- dim(virco)[1]
PrMutSub <- data.frame(PrMut[ , apply(PrMut,2,sum) > NObs*.05])
Trait <- IDV.Fold - NFV.Fold
TtestP <- function(Geno){
	return(t.test(Trait[Geno==1],
		Trait[Geno==0], na.rm=T)$"p.value")
	}
Pvec <- apply(PrMutSub, 2, TtestP)
sort(Pvec)
names(PrMutSub)[Pvec < 0.05]

PvecAdj <- p.adjust(Pvec, method="bonferroni")
sort(PvecAdj)
names(PrMutSub)[PvecAdj < 0.05]

# Example 4.2
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Example 4.2 (Tukeys single-step method):
attach(fms)
Trait <- NDRM.CH
summary(lm(Trait~resistin_c180g))
TukeyHSD(aov(Trait~resistin_c180g))

# Example 4.3
# Reading in Virco data:
virco <- 
read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv",header=T,sep=",")

# Necessary code from Example 4.1:
attach(virco)
PrMut <- virco[,23:121]!="-" & virco[,23:121]!="."
NObs <- dim(virco)[1]
PrMutSub <-data.frame(PrMut[ , apply(PrMut,2,sum) > NObs*.05])
Trait <- IDV.Fold - NFV.Fold
TtestP <- function(Geno){
	return(t.test(Trait[Geno==1],
		Trait[Geno==0], na.rm=T)$"p.value")
	}
Pvec <- apply(PrMutSub, 2, TtestP)
	
# Example 4.3 (Benjamini and Hochberg Adjustment):
Pvec <- as.vector(Pvec)
m <- length(Pvec)
BHp <- sort(Pvec,decreasing=T)*m/seq(m,1)
sort(cummin(BHp))
BHp[order(Pvec,decreasing=T)] <- cummin(BHp)
names(PrMutSub)[BHp < 0.05]
sort(p.adjust(Pvec, method="BH"))

# Example 4.4
# Reading in Virco data:
virco <- 
read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv",header=T,sep=",")

# Necessary code from Example 4.1:
attach(virco)
PrMut <- virco[,23:121]!="-" & virco[,23:121]!="."
NObs <- dim(virco)[1]
PrMutSub <-data.frame(PrMut[ , apply(PrMut,2,sum) > NObs*.05])
Trait <- IDV.Fold - NFV.Fold
TtestP <- function(Geno){
	return(t.test(Trait[Geno==1],
		Trait[Geno==0], na.rm=T)$"p.value")
	}
Pvec <- apply(PrMutSub, 2, TtestP)
Pvec <- as.vector(Pvec)

# Example 4.4 (Benjamini and Yekutieli adjustment)
BYp <- p.adjust(Pvec, method="BY")
sort(BYp)
names(PrMutSub)[BYp < 0.05]

# Example 4.5
# Reading in Virco data:
virco <- 
read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", 
header=T, sep=",")

# Necessary code from Example 4.1:
attach(virco)
PrMut <- virco[,23:121]!="-" & virco[,23:121]!="."
NObs <- dim(virco)[1]
PrMutSub <-data.frame(PrMut[ , apply(PrMut,2,sum) > NObs*.05])
Trait <- IDV.Fold - NFV.Fold
TtestP <- function(Geno){
	return(t.test(Trait[Geno==1],
		Trait[Geno==0], na.rm=T)$"p.value")
	}
Pvec <- apply(PrMutSub, 2, TtestP)
Pvec <- as.vector(Pvec)

# Example 4.5 (Calculation of the q-value)
library(qvalue)
sort(qvalue(Pvec, lambda=0)$qvalues)
sort(qvalue(Pvec, pi0.method="bootstrap")$qvalues)
qvalue(Pvec, pi0.method="bootstrap")$pi0

# Example 4.6
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Example 4.6: Free step-down resampling adjustment
attach(fms)
Actn3Bin <- data.frame(actn3_r577x!="TT", actn3_rs540874!="AA",
				actn3_rs1815739!="TT", 
actn3_1671064!="GG")
Mod <- summary(lm(NDRM.CH~.,data=Actn3Bin))
Mod
TestStatObs <- Mod$coefficients[-1,3]
Tobs <- as.vector(sort(abs(TestStatObs)))
MissDat <- apply(is.na(Actn3Bin),1,any) | is.na(NDRM.CH)
Actn3BinC <- Actn3Bin[!MissDat,]
Ord <- order(abs(TestStatObs))
M <- 1000
NSnps <- 4
Nobs <- sum(!MissDat)
TestStatResamp <- matrix(nrow=M,ncol=NSnps)
for (i in 1:M){ 
		Ynew <- sample(Mod$residuals,size=Nobs,replace=T)
		ModResamp <- summary(lm(Ynew~.,data=Actn3BinC))
		TestStatResamp[i,] <- 
abs(ModResamp$coefficients[-1,3])[Ord]
		}
Qmat <- t(apply(TestStatResamp, 1, cummax))
Padj <- apply(t(matrix(rep(Tobs,M), NSnps)) < Qmat, 2, mean)

# Bootstrap example

set.seed(1)
y1 <- rnorm(75)
mean(y1)-qt(.975, df=74)*sqrt(var(y1)/75)
mean(y1)+qt(.975, df=74)*sqrt(var(y1)/75)
sim <- rep(NA, 1000)
for(i in 1:1000){
  y2 <- sample(y1, replace=TRUE)
  sim[i] <- mean(y2) 
}
quantile(sim, c(.025, .975))
sim <- rep(NA, 10000)
for(i in 1:10000){
  y2 <- sample(y1, replace=TRUE)
  sim[i] <- mean(y2)
}
quantile(sim, c(.025, .975))
 

# Example 4.7
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Neccessary code from Example 4.6:
attach(fms)
Actn3Bin <- data.frame(actn3_r577x!="TT",actn3_rs540874!="AA",
				actn3_rs1815739!="TT",actn3_1671064!="GG")
Mod <- summary(lm(NDRM.CH~.,data=Actn3Bin))
MissDat <- apply(is.na(Actn3Bin),1,sum)>0 | is.na(NDRM.CH)
Actn3BinC <- Actn3Bin[!MissDat,]
Nobs <- sum(!MissDat)

# Example 4.7 (Null unrestricted bootstrap approach):
CoefObs <- as.vector(Mod$coefficients[-1,1])
B <- 1000
NSnps <- 4
Nobs <- sum(!MissDat)
TestStatBoot <- matrix(nrow=B,ncol=NSnps)
for (i in 1:B){ 
		SampID <- sample(1:Nobs,size=Nobs, replace=T)
		Ynew <- NDRM.CH[!MissDat][SampID]
		Xnew <- Actn3BinC[SampID,]
		CoefBoot <- 
summary(lm(Ynew~.,data=Xnew))$coefficients[-1,1]
		SEBoot <- summary(lm(Ynew~.,data=Xnew))$coefficients[-1,2]
		if (length(CoefBoot)==length(CoefObs)){
			TestStatBoot[i,] <- (CoefBoot-CoefObs)/SEBoot
			}
		}
for (cj in seq(2.7,2.8,.01)){
	print(cj)
	print(mean(apply(abs(TestStatBoot)>cj,1,sum)>=1,na.rm=T))
	}
	
# Effective number of tests

corActn3 <- cor(Actn3BinC)
eigenValActn3 <- eigen(corActn3)$values
mEff <- 1+(4-1)*(1-var(eigenValActn3)/4)
mEff
0.05/4
0.05/mEff








