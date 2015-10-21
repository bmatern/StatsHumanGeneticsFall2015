# Ben Matern
# Homework 5

# 1) Foulkes 5.1

fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt",  header=T, sep="\t")
attach(fms)
#install.packages("haplo.stats")
columnNames <- colnames(fms)
sort(columnNames)

library(haplo.stats)

Geno <- cbind(substr(resistin_a537c,1,1), substr(resistin_a537c,2,2),
              substr(resistin_c180g,1,1), substr(resistin_c180g,2,2),
              substr(resistin_c30t,1,1), substr(resistin_c30t,2,2),
              substr(resistin_c398t,1,1), substr(resistin_c398t,2,2),
              substr(resistin_c980g,1,1), substr(resistin_c980g,2,2),
              substr(resistin_g540a,1,1), substr(resistin_g540a,2,2))

Geno
SNPnames <- c("resistin_a537c", "resistin_c180g", "resistin_c30t", "resistin_c398t" , "resistin_c980g",   "resistin_g540a")
SNPnames
Geno.C <- Geno[Race=="Caucasian" & !is.na(Race),]
HaploEM <- haplo.em(Geno.C, locus.label=SNPnames, 
                    control=haplo.em.control(min.posterior=1e-4))
HaploEM

Geno.AA <- Geno[Race=="African Am" & !is.na(Race),]
HaploEM2 <- haplo.em(Geno.AA, locus.label=SNPnames, 
                     control=haplo.em.control(min.posterior=1e-4))
HaploEM2
detach(fms)


#Perform a test to determine if their alleles are different
#See example 5.3
#
# FreqDiff <- HaploEM2$hap.prob[4] - HaploEM$hap.prob[4]
# s1 <- HapFreqSE(HaploEM)[4]
# s2 <- HapFreqSE(HaploEM2)[4]
# SE <- sqrt(s1^2 + s2^2)
# CI <- c(FreqDiff - 1.96*SE, FreqDiff + 1.96*SE)
# CI



# 2) Foulkes 5.2

# Based on the HGDP data, estimate the AKT1 haplotype 
# frequencies within groups defined by the variable Population. 
# Repeat estimation within groups defined by Geographic.origin. 
# Compare and contrast your conclusions for the two analyses.

hgdp <- read.delim("http://people.umass.edu/foulkes/asg/data/HGDP_AKT1.txt",  header=T, sep="\t")
hgdp
library(haplo.stats)
columnNames <- colnames(hgdp)
sort(columnNames)
attach(hgdp)
popNames<-data.matrix(sort(unique(hgdp$Population)))
popNames
geoNames<-data.matrix(sort(unique(hgdp$Geographic.origin)))
geoNames

Geno <- cbind(substr(hgdp$AKT1.C0756A,1,1), substr(hgdp$AKT1.C0756A,2,2),
              substr(hgdp$AKT1.C6024T,1,1), substr(hgdp$AKT1.C6024T,2,2),
              substr(hgdp$AKT1.G2347T,1,1), substr(hgdp$AKT1.G2347T,2,2),
              substr(hgdp$AKT1.G2375A,1,1), substr(hgdp$AKT1.G2375A,2,2))
SNPnames <- c("AKT1.C0756A", "AKT1.C6024T", "AKT1.G2347T", "AKT1.G2375A")


for(i in 1:length(popNames))
{
  populationName<-popNames[i,]
  print(paste(populationName,":"))
  currentPopGeno <- Geno[Population==populationName & !is.na(Population),]
  HaploEM <- haplo.em(currentPopGeno, locus.label=SNPnames, control=haplo.em.control(min.posterior=1e-4))
  print(HaploEM)
}


for(i in 1:length(geoNames))
{
  geoName<-geoNames[i,]
  print(paste(geoName,":"))
  currentGeoGeno <- Geno[Geographic.origin==geoName & !is.na(Geographic.origin),]
  HaploEM <- haplo.em(currentGeoGeno, locus.label=SNPnames, control=haplo.em.control(min.posterior=1e-4))
  print(HaploEM)
}


# 3) Foulkes 5.3
# Is there an association between Gender and Geographic.origin in the HGDP data? 
# If so, how would this influence your interpretation of an analysis of a genotype–trait association?

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

hgdp <- read.delim("http://people.umass.edu/foulkes/asg/data/HGDP_AKT1.txt",  header=T, sep="\t")
hgdp
library(haplo.stats)
columnNames <- colnames(hgdp)
sort(columnNames)
attach(hgdp)


fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", header=T, sep="\t")

# Necessary code from Example 5.1 
#attach(fms)
detach(fms)
#library(haplo.stats)
#Geno <- cbind(substr(actn3_r577x,1,1), substr(actn3_r577x,2,2),
##              substr(actn3_rs540874,1,1), substr(actn3_rs540874,2,2),
#              substr(actn3_rs1815739,1,1), substr(actn3_rs1815739,2,2),
#              substr(actn3_1671064,1,1), substr(actn3_1671064,2,2))
Geographic.origin
Geno <- cbind(substr(Geographic.origin,1,200))
Geno
SNPnames <- c("Geographic.origin")
#Geno.C <- Geno[Race=="Caucasian" & !is.na(Race),]
table(Gender)
Geno.C <- Geno[Gender=="F" & !is.na(Gender),]
Geno.C
HaploEM <- haplo.em(Geno.C,locus.label=SNPnames, 
                    control=haplo.em.control(min.posterior=1e-4))
HaploEM

HapMat <- HapDesign(HaploEM) 
Trait <- NDRM.CH[Race=="Caucasian" & !is.na(Race)]
mod1 <- (lm(Trait~HapMat))
mod2 <- (lm(Trait~1))
anova(mod2,mod1)


# 4) Foulkes 5.4
# Apply haploytpe trend regression (HTR) to determine if there is an association 
# between the resistin haplotypes and change in non-dominant arm muscle strength 
# within African Americans using the FAMuSS data.


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


fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt",  header=T, sep="\t")

attach(fms)
library(haplo.stats)
#Geno <- cbind(substr(actn3_r577x,1,1), substr(actn3_r577x,2,2),
#              substr(actn3_rs540874,1,1), substr(actn3_rs540874,2,2),
#              substr(actn3_rs1815739,1,1), substr(actn3_rs1815739,2,2),
#              substr(actn3_1671064,1,1), substr(actn3_1671064,2,2))
Geno <- cbind(substr(resistin_a537c,1,1), substr(resistin_a537c,2,2),
              substr(resistin_c180g,1,1), substr(resistin_c180g,2,2),
              substr(resistin_c30t,1,1), substr(resistin_c30t,2,2),
              substr(resistin_c398t,1,1), substr(resistin_c398t,2,2),
              substr(resistin_c980g,1,1), substr(resistin_c980g,2,2),
              substr(resistin_g540a,1,1), substr(resistin_g540a,2,2))
#SNPnames <- c("actn3_r577x", "actn3_rs540874", "actn3_rs1815739",
#              "actn3_1671064")
SNPnames <- c("resistin_a537c", "resistin_c180g", "resistin_c30t", "resistin_c398t" , "resistin_c980g",   "resistin_g540a")
Geno.C <- Geno[Race=="African Am" & !is.na(Race),]
#table(Race)
HaploEM <- haplo.em(Geno.C,locus.label=SNPnames, control=haplo.em.control(min.posterior=1e-4))

HapMat <- HapDesign(HaploEM) 
Trait <- NDRM.CH[Race=="African Am" & !is.na(Race)]
mod1 <- (lm(Trait~HapMat))
mod2 <- (lm(Trait~1))
anova(mod2,mod1)

# 5) Foulkes 5.5
# Using the expectation-maximization approach of the haplo.glm() func- tion, 
# determine if there is an association between the resistin haplotypes 
# and change in non-dominant arm muscle strength, as measured by NDRM.CH
# , within African Americans, based on the FAMuSS data. 
# Consider both domi- nant and additive genetic models.


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


fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt",  header=T, sep="\t")


Geno <- cbind(substr(resistin_a537c,1,1), substr(resistin_a537c,2,2),
              substr(resistin_c180g,1,1), substr(resistin_c180g,2,2),
              substr(resistin_c30t,1,1), substr(resistin_c30t,2,2),
              substr(resistin_c398t,1,1), substr(resistin_c398t,2,2),
              substr(resistin_c980g,1,1), substr(resistin_c980g,2,2),
              substr(resistin_g540a,1,1), substr(resistin_g540a,2,2))
#SNPnames <- c("actn3_r577x", "actn3_rs540874", "actn3_rs1815739",
#              "actn3_1671064")
SNPnames <- c("resistin_a537c", "resistin_c180g", "resistin_c30t", "resistin_c398t" , "resistin_c980g",   "resistin_g540a")
Geno.C <- Geno[Race=="African Am" & !is.na(Race),]


library(haplo.stats)
Geno.C <- setupGeno(Geno.C)
Trait <- NDRM.CH[Race=="African Am" & !is.na(Race)]
Dat <- data.frame(Geno.C=Geno.C, Trait=Trait)
haplo.glm(Trait~Geno.C,data=Dat, 
          allele.lev=attributes(Geno.C)$unique.alleles)
haplo.glm(Trait~Geno.C,data=Dat, 
          allele.lev=attributes(Geno.C)$unique.alleles, 
          control=haplo.glm.control(haplo.base=9))
haplo.glm(Trait~Geno.C,data=Dat, 
          allele.lev=attributes(Geno.C)$unique.alleles,
          control=haplo.glm.control(haplo.effect="dominant"))


