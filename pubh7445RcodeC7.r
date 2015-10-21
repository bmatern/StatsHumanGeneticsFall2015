# Example 7.1
# Reading in Virco data:
virco <- 
read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", 
header=T, sep=",")

# Example 7.1 (An application of random forests):
install.packages("randomForest")
library(randomForest)
attach(virco)
Trait <- NFV.Fold - IDV.Fold
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
Trait.c <- Trait[!is.na(Trait)]
VircoGeno.c <- VircoGeno[!is.na(Trait),]
RegRF <- randomForest(VircoGeno.c, Trait.c, importance=TRUE)
RegRF
varImpPlot(RegRF,main="")

# Example 7.2
# Reading in FAMuSS data:
fms <- 
read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt",header=T,sep="\t")

# Example 7.2 (RF with missing SNP data - single imputation):
attach(fms)
Trait <- NDRM.CH[Race=="Caucasian" & !is.na(Race) & !is.na(NDRM.CH)]
NamesAkt1Snps <- names(fms)[substr(names(fms),1,4)=="akt1"]
FMSgeno <- fms[,is.element(names(fms),NamesAkt1Snps)][Race=="Caucasian" & 
!is.na(Race) &!is.na(NDRM.CH),]
dim(FMSgeno)
round(apply(is.na(FMSgeno),2,sum)/dim(FMSgeno)[1],3)
library(randomForest)
FMSgenoRough <- na.roughfix(FMSgeno)
table(FMSgeno$"akt1_t22932c")
round(apply(is.na(FMSgenoRough),2,sum)/dim(FMSgeno)[1],3)
RandForRough <- randomForest(FMSgenoRough,Trait,importance=TRUE)
RandForRough$"importance"[order(RandForRough$"importance"[,1],decreasing=TRUE),]

# Example 7.3
# Reading in FAMuSS data:
fms <- 
read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt",header=T,sep="\t")

# Necessary code from Example 7.2
attach(fms)
Trait <- NDRM.CH[Race=="Caucasian" & !is.na(Race) & !is.na(NDRM.CH)]
NamesAkt1Snps <- names(fms)[substr(names(fms),1,4)=="akt1"]
FMSgeno <- fms[,is.element(names(fms),NamesAkt1Snps)][Race=="Caucasian" & 
!is.na(Race) &!is.na(NDRM.CH),]
library(randomForest)
FMSgenoRough <- na.roughfix(FMSgeno)

# Example 7.3 (RF with missing SNP data - multiple imputation):
Trait <- NDRM.CH[Race=="Caucasian" & !is.na(Race) & !is.na(NDRM.CH)]
NamesAkt1Snps <- names(fms)[substr(names(fms),1,4)=="akt1"]
FMSgeno <- fms[,is.element(names(fms), NamesAkt1Snps)][Race=="Caucasian" & 
!is.na(Race) & !is.na(NDRM.CH),]
FMSgenoMI <- rfImpute(FMSgeno, Trait)
RandForFinal <- randomForest(FMSgenoMI[,-1], Trait, importance=TRUE)
RandForFinal$"importance"[order(RandForFinal$"importance"[,1], 
decreasing=TRUE),]
table(FMSgenoMI$akt1_t10726c_t12868c, FMSgenoRough$akt1_t10726c_t12868c)

# Example 7.4

The package mirf is no longer supported

# Example 7.5
# Reading in Virco data:
virco <- 
read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", 
header=T, sep=",")

# Example 7.5 (Application of logic regression):
attach(virco)
Trait <- NFV.Fold - IDV.Fold
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
Trait.c <- Trait[!is.na(Trait)]
VircoGeno.c <- VircoGeno[!is.na(Trait),]
install.packages("LogicReg")
library(LogicReg)
VircoLogicReg <- logreg(resp=Trait.c, bin=VircoGeno.c, select=1)
plot(VircoLogicReg)
VircoLogicReg
VircoLogicRegMult <- logreg(resp=Trait.c, bin=VircoGeno.c, select=2, 
ntrees=2, nleaves=8)
plot(VircoLogicRegMult)
VircoLogicRegMult

# Example 7.6
# Reading in Virco data:
virco <- 
read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", 
header=T, sep=",")

# Example 7.6 (Monte Carolo logic regression):
library(LogicReg)
attach(virco)
Trait <- SQV.Fold
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
Trait.c <- Trait[!is.na(Trait)]
VircoGeno.c <- VircoGeno[!is.na(Trait),]
VircoLogicRegMCMC <- logreg(resp=Trait.c, bin=VircoGeno.c, select=7)
plot(sort(VircoLogicRegMCMC$single), xlab="Sorted SNPs", ylab="Number of 
selected models")
names(VircoGeno)[order(VircoLogicRegMCMC$single)]

# Example 7.7
# Reading in Virco data:
virco <- 
read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", 
header=T, sep=",")

# Example 7.7 (An application of MARS):
attach(virco)
Trait <- NFV.Fold - IDV.Fold
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
Trait.c <- Trait[!is.na(Trait)]
VircoGeno.c <- VircoGeno[!is.na(Trait),]
install.packages("earth")
library(earth)
VircoMARS <- earth(Trait.c~., data=VircoGeno.c, degree=2)
summary(VircoMARS)
evimp(VircoMARS)

