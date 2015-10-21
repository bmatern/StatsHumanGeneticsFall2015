# Ben Matern
# Homework 4

# 1)
# T test
pval <- rep(NA,1000)
for(i in 1:1000)
{
  y1 <- rnorm(40, mean=1, sd=2)
  y2 <- rnorm(40, mean=2, sd=2)
  pval[i] <- t.test(y1, y2, var.equal=TRUE)$p.value
}
sum(pval<.05)/1000

# Wilcoxan 
pval <- rep(NA,1000)
for(i in 1:1000)
{
  y1 <- rnorm(40, mean=1, sd=2)
  y2 <- rnorm(40, mean=2, sd=2)
  pval[i] <- wilcox.test(y1, y2, var.equal=TRUE)$p.value
}
sum(pval<.05)/1000

# 2)

# Foulkes 4.2
#virco <- read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", header=T, sep=",")
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Example 4.1 (Bonferroni adjustment):
#attach(virco)
#attach(fms)
#columns <- colnames(fms)
#sortCols <- sort(columns)
#columns
#sortCols
#HDL_C = 312
#AKT = 15:37
#virco
fms
#PrMut <- virco[,23:121]!="-" & virco[,23:121]!="."
submatrix <- fms[,15:37]
PrMut <- submatrix
#calculate the minor allele and see if each column has one.
for(i in 1:23)
{
  GenoCount <- table(submatrix[i])
  NumbObs <- sum(!is.na(submatrix[i]))
  if(GenoCount[1] > GenoCount[3] || is.na(GenoCount[3]))
  {
    #first allele is the major allele, second is minor/variant
    minorAllele <- substr(names(GenoCount)[3],1,1)
  }
  else
  {
    minorAllele <- substr(names(GenoCount)[1],1,1)
  }
  
  if(!is.na(minorAllele))
  {
    PrMut[,i] <- (substr(submatrix[,i],1,1) == minorAllele)
  }
  else
  {
    #hello<-"hello"
    #hello
    #nothing has the minor allele.
    PrMut[,i] <- (substr(submatrix[,i],1,1) != substr(submatrix[,i],2,2))
  }
}
#NA doesnt have th allele.
PrMut[is.na(PrMut)] <- FALSE
# PrMut[is.na(PrMut)] 


PrMut
NObs <- dim(fms)[1]
PrMutSub <- data.frame(PrMut[ , apply(PrMut,2,sum) > NObs*.05])
PrMutSub
#Trait <- IDV.Fold - NFV.Fold
Trait <- fms$HDL_C
Trait
TtestP <- function(Geno){return(t.test(Trait[Geno==1], Trait[Geno==0], na.rm=T)$"p.value")}



Pvec <- apply(PrMutSub, 2, TtestP)
sort(Pvec)
names(PrMutSub)[Pvec < 0.05]

PvecAdj <- p.adjust(Pvec, method="bonferroni")
sort(PvecAdj)
names(PrMutSub)[PvecAdj < 0.05]
#detach(virco)




Pvec <- apply(PrMutSub, 2, TtestP)
Pvec <- as.vector(Pvec)
m <- length(Pvec)
BHp <- sort(Pvec,decreasing=T)*m/seq(m,1)
sort(cummin(BHp))
BHp[order(Pvec,decreasing=T)] <- cummin(BHp)
names(PrMutSub)[BHp < 0.05]
sort(p.adjust(Pvec, method="BH"))
# detach(virco)




#detach(fms)




# Necessary code from Example 4.1:
attach(virco)
PrMut <- virco[,23:121]!="-" & virco[,23:121]!="."
NObs <- dim(virco)[1]
PrMutSub <-data.frame(PrMut[ , apply(PrMut,2,sum) > NObs*.05])
Trait <- IDV.Fold - NFV.Fold
Trait
TtestP <- function(Geno){return(t.test(Trait[Geno==1],Trait[Geno==0], na.rm=T)$"p.value")}
Pvec <- apply(PrMutSub, 2, TtestP)

# Example 4.3 (Benjamini and Hochberg Adjustment):
Pvec <- as.vector(Pvec)
m <- length(Pvec)
BHp <- sort(Pvec,decreasing=T)*m/seq(m,1)
sort(cummin(BHp))
BHp[order(Pvec,decreasing=T)] <- cummin(BHp)
names(PrMutSub)[BHp < 0.05]
sort(p.adjust(Pvec, method="BH"))
detach(virco)




# 3)

# Foulkes 4.3
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")
#attach(fms)
#Trait <- NDRM.CH
Trait <- fms$CHOL
Trait
resistin_c180g
#lm is linear model
summary(lm(Trait~fms$resistin_c180g))
TukeyHSD(aov(Trait~fms$resistin_c180g))


#colnames(fms)
#230 = gender
#193 = resistin_c180g
#CHOL = 311
#Strat by Gender
noGender<-fms[is.na(fms$Gender),c(230,193,311)]
#trim the nulls
trimmeddata <- fms[!is.na(fms$Gender),]


maleData<-trimmeddata[trimmeddata$Gender=="Male",c(230,193,311)]
#maleData<-maleData[!is.na(maleData$Gender),]
maleTrait<-maleData$CHOL
summary(lm(maleTrait~maleData$resistin_c180g))
TukeyHSD(aov(maleTrait~maleData$resistin_c180g))

femaleData<-trimmeddata[trimmeddata$Gender=="Female",c(230,193,311)]
femaleTrait<-femaleData$CHOL
summary(lm(femaleTrait~femaleData$resistin_c180g))
TukeyHSD(aov(femaleTrait~femaleData$resistin_c180g))


#detach(fms)

# 4 )
# foulkes 4.6
# Example 4.6: Free step-down resampling adjustment;
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")
#attach(fms)
sort(colnames(fms))

#Find upper quartile
trimmedData<-fms[!is.na(fms$NDRM.CH),]
trimmedData$NDRM.CH
ndrmdata<-sort(trimmedData$NDRM.CH)
round(length(ndrmdata)*3/4)
qCutoff<-ndrmdata[round(length(ndrmdata)*3/4)]
isUpperQuartile<-trimmedData$NDRM.CH > qCutoff
isUpperQuartile

Actn3Bin <- data.frame(trimmedData$actn3_r577x!="TT", trimmedData$actn3_rs540874!="AA",
                       trimmedData$actn3_rs1815739!="TT", 
                       trimmedData$actn3_1671064!="GG")

Mod <- summary(lm(isUpperQuartile~.,data=Actn3Bin))
Mod
TestStatObs <- Mod$coefficients[-1,3]
TestStatObs
Tobs <- as.vector(sort(abs(TestStatObs)))
Tobs
MissDat <- apply(is.na(Actn3Bin),1,any) | is.na(trimmedData$NDRM.CH)
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
Padj

# 5 ) 


# Effective number of tests
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")
#attach(fms)
#Actn3Bin <- data.frame(fms$actn3_r577x!="TT", fms$actn3_rs540874!="AA",
#                       fms$actn3_rs1815739!="TT", 
#                      fms$actn3_1671064!="GG")
Mod <- summary(lm(fms$NDRM.CH~.,data=Actn3Bin))
Mod
TestStatObs <- Mod$coefficients[-1,3]
Tobs <- as.vector(sort(abs(TestStatObs)))
MissDat <- apply(is.na(Actn3Bin),1,any) | is.na(fms$NDRM.CH)
Actn3BinC <- Actn3Bin[!MissDat,]

corActn3 <- cor(Actn3BinC)
eigenValActn3 <- eigen(corActn3)$values
mEff <- 1+(4-1)*(1-var(eigenValActn3)/4)
mEff
0.05/4
0.05/mEff

#For real this time 
hgdp <- read.delim("http://people.umass.edu/foulkes/asg/data/HGDP_AKT1.txt",  header=T, sep="\t")
colnames(hgdp)
data<-rnorm(1064)

Actn3Bin <- data.frame(hgdp$AKT1.C0756A!="CC", 
                       hgdp$AKT1.C6024T!="CC",
                       hgdp$AKT1.G2347T!="GG", 
                       hgdp$AKT1.G2375A!="GG")

Mod <- summary(lm(data~.,data=Actn3Bin))
Mod
TestStatObs <- Mod$coefficients[-1,3]
Tobs <- as.vector(sort(abs(TestStatObs)))
MissDat <- apply(is.na(Actn3Bin),1,any) 
Actn3BinC <- Actn3Bin[!MissDat,]

corActn3 <- cor(Actn3BinC)
eigenValActn3 <- eigen(corActn3)$values
mEff <- 1+(4-1)*(1-var(eigenValActn3)/4)
mEff
0.05/4
0.05/mEff



# 6 )

# Independently simulate binary indicator variables for 
# 2000 markers for 100 subjects with success probability 0.5. 
simMatrix<-matrix(NA, 2000, 100)
for(j in 1:100)
{
  for(i in 1:2000)
  {
    #probability  = .5
    simMatrix[i,j]<-(rnorm(1) > 0)
  }
}

#now sim matrix contains a 2000x100 matrix of booleans

outcomes<-rep(0,100)


# for each subject simulate a normally distributed outcome variable that depends on 
# the first 10 markers with regression coefficients of 10,9,...,1 and has standard deviation 1.

# what?  how to i enforce the standard deviation? rnorm?
# Am i sampling from a normal distribution?  
for(j in 1:100)
{
  outcomes[j] = 10 * simMatrix[1,j] + 
    9 * simMatrix[2,j] +
    8 * simMatrix[3,j] +
    7 * simMatrix[4,j] +
    6 * simMatrix[5,j] +
    5 * simMatrix[6,j] +
    4 * simMatrix[7,j] +
    3 * simMatrix[8,j] +
    2 * simMatrix[9,j] +
    1 * simMatrix[10,j]
}

# For each marker, use a 2 sample t-test with equal variance 
# to test if the marker is associated with the trait.
pvalues = rep(0,2000)
for(i in 1:2000)
{
  markerTrait <- simMatrix[i,]
  #g1<-data[whp, ALL_bcrneg$mol.biol=="BCR/ABL"]
  #g2<-data[whp,ALL_bcrneg$mol.biol=="NEG"]
  summt<-t.test(markerTrait, outcomes)
  tobs<-summt$p.value
  pvalues[i]<-tobs
  
}
pvalues
hist(pvalues)

qvalues <- rep(0,2000)
#install.packages("qvalue")
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")


library(qvalue)
pvalues
qvalues <- sort(qvalue(pvalues, lambda=0)$qvalues)
qvalues
sort(qvalue(pvalues, pi0.method="bootstrap")$qvalues)
qvalue(pvalues, pi0.method="bootstrap")$pi0
# 7

set.seed(1)
sim1 <- rnorm(100, mean=0, sd=1)
sim2 <- rnorm(50, mean=3, sd=.5)
combinedSim <- c(sim1, sim2)

hist(combinedSim)

mean(combinedSim)-qt(.975, df=74)*sqrt(var(combinedSim)/75)
mean(combinedSim)+qt(.975, df=74)*sqrt(var(combinedSim)/75)
#sim contains the bootstrap values/
sim <- rep(NA, 1000)
for(i in 1:1000){
  y2 <- sample(combinedSim, replace=TRUE)
  sim[i] <- mean(y2) 
}
quantile(sim, c(.025, .975))
