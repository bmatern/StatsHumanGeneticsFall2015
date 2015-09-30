# Ben Matern
# Homework 2.

# 1)
# 6.5
diceRolls <- sample(1:6, 100, replace=T)
sum(diceRolls==6)

# 6.6
ballLottery <- sample(1:49, 6)
max(ballLottery)
min(ballLottery)

# 6.7
qnorm(.05, mean=0, sd=1)

# 6.8
qnorm(.5 + .05/2, mean=0, sd=1) 

# 6.9
pnorm( 1.5, mean=0, sd=2)

# 2)
# 2.2
setwd("/Users/bmatern/school/Fall2015/Stats4HumGenomics/StatsHumanGeneticsFall2015/Hwk2")
fmsURL<-"http://people.umass.edu/foulkes/asg/data/FMS_data.txt"
fms<-read.delim(file=fmsURL, header=TRUE, sep="\t")
attach(fms)
install.packages("genetics", depend=TRUE)
library(genetics)

NamesAkt2Snps <- names(fms)[substr(names(fms),1,4)=="akt2"]
fmsAkt2 <- fms[,is.element(names(fms),NamesAkt2Snps)]
TtestPval <- function(Geno)
{
  alleleMajor <- allele.names(genotype(Geno, sep="", reorder="freq"))[1]
  GenoWt <- paste(alleleMajor, alleleMajor, sep="")
  GenoBin <- as.numeric(Geno!=GenoWt)[!is.na(Geno)]
  Trait <- NDRM.CH[!is.na(Geno)]
  return(t.test(Trait[GenoBin==1], Trait[GenoBin==0])$p.value)
}

apply(fmsAkt2,2,TtestPval)

# 3)
# 2.4

setwd("/Users/bmatern/school/Fall2015/Stats4HumGenomics/StatsHumanGeneticsFall2015/Hwk2")
fmsURL<-"http://people.umass.edu/foulkes/asg/data/FMS_data.txt"
fms<-read.delim(file=fmsURL, header=TRUE, sep="\t")
attach(fms)
library(genetics)

dim(fms)
str(fms[,1:10])
#Geno<-esr1_rs1042717
Geno<-resistin_c180g 
Geno
trait<-as.numeric(pre.BMI > 25)
plot(Geno)

obsTab <- table(trait,Geno)
obsTab

homozygousVariantOdds <- obsTab[2,3]/obsTab[1,3]
homozygousWildtypeOdds <- obsTab[2,1]/obsTab[1,1]

variantToWildtype <- homozygousVariantOdds / homozygousWildtypeOdds
variantToWildtype

# 4)
fmsURL<-"http://people.umass.edu/foulkes/asg/data/FMS_data.txt"
fms<-read.delim(file=fmsURL, header=TRUE, sep="\t")
attach(fms)
Geno <- actn3_r577x
Trait <- NDRM.CH
ModFull <- lm(Trait~Geno+Gender+Geno*Gender, na.action=na.exclude)
summary(ModFull)
lm(formula=Trait ~ Geno + Gender + Geno * Gender, na.action=na.exclude)
ModReduced <- lm(Trait~Geno+Gender, na.action=na.exclude)
anova(ModReduced, ModFull)

# 5) 

numObs <- 1000
numAA <- 649
numAB <- 300
numBB <- 51
freqA <- (numAA*2 + numAB) / (2*numObs)
freqB <- (numAB + numBB*2) / (2*numObs)

obsCount <- c(numAA, numAB, numBB)
expCount <- c(numObs*freqA^2, 2*numObs*freqA*(freqB), numObs*(freqB)^2)

ChiSqStat <- sum((obsCount - expCount)^2/expCount)
ChiSqStat

numObs <- 1000
numAA <- 640
numAB <- 360
numBB <- 0
freqA <- (numAA*2 + numAB) / (2*numObs)
freqB <- (numAB + numBB*2) / (2*numObs)

obsCount <- c(numAA, numAB, numBB)
expCount <- c(numObs*freqA^2, 2*numObs*freqA*(freqB), numObs*(freqB)^2)

ChiSqStat <- sum((obsCount - expCount)^2/expCount)
ChiSqStat

# 6) 
# a)
numSamples = 200
sampleMeans <- rep(0,numSamples)
for(i in 1:numSamples)
{
  samples <- rnorm(numSamples, 3, 1.5)
  sampleMeans[i] <- mean(samples)
}
plot(sampleMeans)
isGreaterPercentage = 100 *sum(sampleMeans >= 3.2) / numSamples

# b)
numSamples = 100
sampleMeans <- rep(0,numSamples)
for(i in 1:numSamples)
{
  samples <- rnorm(numSamples, 3, 1.5)
  sampleMeans[i] <- mean(samples)
}
plot(sampleMeans)
isGreaterPercentage = 100 *sum(sampleMeans >= 3.2) / numSamples

# 7)
strainA <- c(132,72, 102, 115, 59, 103, 86, 159, 60, 94, 80, 97)
meanA <- mean(strainA)
sdA <- sd(strainA)
strainB <- c(101, 96, 93, 106, 81, 77, 106, 97, 74)
meanB <- mean(strainB)
sdB <- sd(strainB)

# a)
testResults = t.test(strainA, strainB)
testResults

# b) 
# c) 
wcxResults = wilcox.test(strainA, strainB)
wcxResults

# d)
perm <- 1000
tstar <- rep(NA,perm)
for(i in 1:perm)
{
  g1 <- sample(strainA,12)
  g2 <- sample(strainB,9)
  tstar[i] = t.test(g1,g2)$statistic
}
tstar
plot(density(tstar), main= "permuted distribution")
pvalue <- mean(abs(tstar) >= abs(perm))
pvalue

#8 )

#9 )

setwd("/Users/bmatern/school/Fall2015/Stats4HumGenomics/StatsHumanGeneticsFall2015/Hwk2")
fmsURL<-"http://people.umass.edu/foulkes/asg/data/FMS_data.txt"
fms<-read.delim(file=fmsURL, header=TRUE, sep="\t")
attach(fms)
library(genetics)
Geno <- as.factor(actn3_r577x)
# Geno <- as.factor(resistin_c180g)
Trait <- NDRM.CH

# a)
stripchart(Trait ~ Geno)

# b)
AnovaMod <- lm(Trait~Geno, na.action=na.exclude)
AnovaMod

# c)
maxTrait <- max(Trait, na.rm=TRUE)
indexMax <- which(Trait==maxTrait)
modifiedTrait <- Trait[-indexMax]
modifiedGeno <- Geno[-indexMax]
#modifiedTrait <- Trait[! Trait %in% maxTrait]
maxTrait <- max(modifiedTrait, na.rm=TRUE)
indexMax <- which(modifiedTrait==maxTrait)
modifiedTrait <- modifiedTrait[-indexMax]
modifiedGeno <- modifiedGeno[-indexMax]

stripchart(modifiedTrait ~ modifiedGeno)
AnovaMod <- lm(modifiedTrait~modifiedGeno, na.action=na.exclude)
AnovaMod

# d)
kruskal.test(Trait, Geno, na.action=na.exclude)








