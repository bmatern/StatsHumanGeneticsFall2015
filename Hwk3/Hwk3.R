# Hwk3
# Ben Matern

# 1)
# Foulkes 1.3)
# actn3 1671064 gene
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

GenoCount <- table(actn3_1671064)
GenoCount
NumbObs <- sum(!is.na(actn3_1671064))
NumbObs
GenoFreq <- as.vector(GenoCount/NumbObs)
GenoFreq
FreqA <- (2*GenoFreq[1] + GenoFreq[2])/2
FreqA
FreqG <- (GenoFreq[2] + 2*GenoFreq[3])/2
FreqG

#Stratify by Race
TableRaces <- table(Race)
TableRaces

GenoCount <- table(actn3_1671064, Race)
GenoCount
NumbObs <- sum(!is.na(actn3_1671064))
NumbObs
columnSums <- colSums(GenoCount)
columnSums
GenoFreq <- GenoCount
Afrequencies <- rep(0,6)
Gfrequencies <- rep(0,6)
for(i in 1:6)
{
  for (j in 1:3 ) 
  {
     GenoFreq[j + (i-1) * 3] <- GenoCount[j + (i-1) * 3] / columnSums[i]
     
  }
  Afrequencies[i] <- (GenoCount[(i-1) * 3 + 1] * 2 + GenoCount[(i-1) * 3 + 2]) / (columnSums[i] * 2)
  Gfrequencies[i] <- 1 - Afrequencies[i]
}

columnSums
Afrequencies
Gfrequencies
GenoCount
GenoFreq

#2. Foulkes 1.5

vircoURL <- "http://people.umass.edu/foulkes/asg/data/Virco_data.csv"
virco <- read.csv(file=vircoURL)
dim(virco)
#virco[1:5,c(1,6,11,32,85,93,104,112,122)]
#virco[1:5,104]
mutationData <- virco[,c(23,32,52,93,104,112)]
mutationData
dim(mutationData)
totalCount <- rep(dim(mutationData)[1],6)

#numberWildType <- sum(virco[])
columnSums <- colSums(mutationData=="-")
wildTypeProportions <- columnSums / totalCount
mutationProportions <- 1 - wildTypeProportions 

mutationProportions

#P1, P10, P30, P71, P82 and P90.

#3 Foulkes 3.2
hgdp <-read.delim("http://people.umass.edu/foulkes/asg/data/HGDP_AKT1.txt",  header=T, sep="\t")
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")
#detachAllData()
hgdp
fms
head(fms)
head(hgdp)

attach(fms)
Akt1Snp1 <- akt1_t10726c_t12868c
#akt1 t10726c t12868c
#akt1_t10726c_t12868c
ObsCount <- table(Akt1Snp1)


# More than 5 observations, do the chi square method.

Nobs <- sum(ObsCount)
ObsCount
FreqC <- (2 * ObsCount[1] + ObsCount[2])/(2*Nobs)
FreqC
ExpCount <- c(Nobs*(FreqC)^2, 2*Nobs*FreqC*(1-FreqC),Nobs*(1-FreqC)^2)
ExpCount
ChiSqStat <- sum((ObsCount - ExpCount)^2/ExpCount)
ChiSqStat
qchisq(1-0.05,df=1)

#Stratify by Race
Strat <- table(akt1_t10726c_t12868c, Race)
Nobs <- sum(Strat)
head(Strat)
columnSums <- colSums(Strat)

FreqCStrat <- rep(0,6)
ChiSqStatStrat <- rep(0,6)
for(i in 1:6)
{  
  FreqCStrat[i] <- (Strat[1 + (i-1)*3] * 2 + Strat[2 + (i-1)*3])/(2*columnSums[i])
  ObsCount <- c(Strat[1+ (i-1) * 3], Strat[2+ (i-1) * 3], Strat[3+ (i-1) * 3])
  FreqC <- FreqCStrat[i]
  ExpCount <- c(columnSums[i]*(FreqC)^2, 2*columnSums[i]*FreqC*(1-FreqC),columnSums[i]*(1-FreqC)^2)
  ObsCount
  ExpCount
  #Strat[1]*2 + 
  ChiSqStatStrat[i] <- sum((ObsCount - ExpCount)^2/ExpCount)
}
FreqCStrat
ChiSqStatStrat


#4 Foulkes 3.4
#Report estimates of pairwise linkage disequilibrium (LD) for all SNPs within the akt1 gene 
#for the HGDP data. Do these estimates tend to vary across Geographic.area? 
#Interpret your findings.

hgdp <-read.delim("http://people.umass.edu/foulkes/asg/data/HGDP_AKT1.txt",  header=T, sep="\t")

fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt",  header=T, sep="\t")

head(fms)
head(hgdp)

attach(hgdp)
library(genetics)


#AKT1.C0756A AKT1.C6024T AKT1.G2347T  AKT1.G2375A
Actn3Snp1 <- genotype(AKT1.C0756A,sep="")
Actn3Snp2 <- genotype(AKT1.C6024T,sep="")
Actn3Snp3 <- genotype(AKT1.G2347T,sep="")
Actn3Snp4 <- genotype(AKT1.G2375A,sep="")
Actn3AllSnps <- data.frame(Actn3Snp1,Actn3Snp2,Actn3Snp3,Actn3Snp4)
LD(Actn3AllSnps)$"D'"
#install.packages("LDheatmap")
library(LDheatmap)
LDheatmap(Actn3AllSnps, LDmeasure="D'")

#5 Foulkes 3.6
#Determine whether there is any evidence for population substructure in 
#African Americans in the FAMuSS data. Explain how you reached this conclusion.
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

#table(Geographic.area)
##head(fms)
#library(genetics)
#Akt1Snp1 <- genotype(AKT1.C0756A, sep="")
#HWEGeoArea <- tapply(Akt1Snp1,INDEX=Geographic.area,HWE.chisq)  
#HWEGeoArea$"Central Africa"
#HWEGeoArea$"South America"
#HWEGeoArea$""

table(Race)
subset = fms[fms$Race=="African Am",]
head(subset)
#dat<-fms[!is.na(fms$pre.BMI) & !is.na(fms$actn3_rs540874) , index]
NamesAkt1Snps <- names(subset)[substr(names(subset),1,4)=="akt1"]
FMSgeno <- subset[,is.element(names(subset),NamesAkt1Snps)]
FMSgenoNum <- data.matrix(FMSgeno)
FMSgenoNum[is.na(FMSgenoNum)] <- 4
DistFmsGeno <- as.matrix(dist(FMSgenoNum))

plot(cmdscale(DistFmsGeno),xlab="C1",ylab="C2")
#abline(v=0,lty=2)	
#abline(h=4,lty=2)
