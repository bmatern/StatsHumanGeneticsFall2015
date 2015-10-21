#####################################################
## setting up data

attach(fms)

akt1Dat <- fms[,substr(names(fms),1,4)=="akt1"]
#> dim(akt1Dat)
#[1] 1397   24

# check if all SNPs have 3 levels

q1=rep(NA,24)
for(i in 1:24) q1[i]=nlevels(akt1Dat[,i])

which(q1!=3)

table(akt1Dat[,15])
table(akt1Dat[,17])

# exclude the 15th SNP and make the Und elements of akt1Dat into NA

levels(akt1Dat[,17])=c("AA", "GA", "GG",NA)

akt1Dat=akt1Dat[,-15]

# here is indicator of having at least 1 copy of the minor allele

getMnrAll <- function(y){
  t1=table(y)
  n1=levels(y)
  if(t1[1]<=t1[3]) x=ifelse(y!=n1[3],1,0)
  else x=ifelse(y!=n1[1],1,0)
  x
}  

akt1Bin=matrix(NA,1397,23)
for(i in 1:23) akt1Bin[,i]=getMnrAll(akt1Dat[,i])

# get a complete data set

MissDat <- apply(is.na(akt1Bin),1,any) | is.na(NDRM.CH)
akt1BinC <- akt1Bin[!MissDat,]

> dim(akt1BinC)
[1] 707  23

TraitC <- NDRM.CH[!MissDat]

# now dichotomize

TraitCD <- TraitC > median(TraitC)

> summary(TraitCD)
   Mode   FALSE    TRUE    NA's 
logical     398     309       0 

then set 100 observations aside as test cases

set.seed(102)
testSmp <- sample(1:707,100)
akt1BinC.ts <- akt1BinC[testSmp,]
Trait.ts <- TraitCD[testSmp]
akt1BinC.ls <- akt1BinC[-testSmp,]
Trait.ls <- TraitCD[-testSmp]

> summary(Trait.ts)
   Mode   FALSE    TRUE    NA's 
logical      56      44       0 
> summary(Trait.ls)
   Mode   FALSE    TRUE    NA's 
logical     342     265       0 

# so 44% TRUE overall and in both samples

######################################
## linear model

clsLM <- lm(Trait.ls~.,data=data.frame(akt1BinC.ls))
prdLM <- predict(clsLM,newdata=data.frame(akt1BinC.ts))

install.packages("pROC")
library(pROC)

roc(Trait.ts,prdLM,plot=T)

#############################
## random forests

install.packages("randomForest")
library(randomForest)

clsRF <- randomForest(akt1BinC.ls, factor(Trait.ls))

prdRF <- predict(clsRF,akt1BinC.ts)

t1=table(prdRF,Trait.ts)
#       Trait.ts
#prdRF   FALSE TRUE
#  FALSE    50   39
#  TRUE      6    5

t1[2,2]/sum(t1[,2])
t1[1,1]/sum(t1[,1])

prdRFp <- predict(clsRF,akt1BinC.ts,type="prob")

install.packages("pROC")
library(pROC)

rfROC <- roc(Trait.ts,prdRFp[,1],plot=TRUE)

# but the output is random, so perhaps average over many tries?

# hard way: by hand

#M=10
#sensRF=matrix(NA,M,100)
#specRF=matrix(NA,M,100)
#for(i in 1:M){
#  clsRF <- randomForest(akt1BinC.ls, factor(Trait.ls))
#  prdRFp <- predict(clsRF,akt1BinC.ts,type="prob")
#  rfROC <- roc(Trait.ts,prdRFp[,1])
#  nn <- length(rfROC$sensitivities)
#  sensRF[i,1:nn]=rfROC$sensitivities
#  specRF[i,1:nn]=rfROC$specificities
#}

#plot(specRF[1,],sensRF[1,],type="l")
#for(i in 2:10) lines(specRF[i,],sens[i,])

# or simpler

for(i in 1:M){
  clsRF <- randomForest(akt1BinC.ls, factor(Trait.ls))
  prdRFp <- predict(clsRF,akt1BinC.ts,type="prob")
  rfROC <- lines.roc(Trait.ts,prdRFp[,1])
}

############################
## MARS

install.packages("earth")
library(earth)   

# a little tricky how predict works with earth, for example

clsMARS <- earth(factor(Trait.ls)~akt1BinC.ls, degree=2)

prdMARS <- predict(clsMARS,akt1BinC.ts)

# gives an error, while the following works, although not sure how to set threshold, default is 0.5 but that performs poorly

clsMARS <- earth(factor(Trait.ls)~., data=data.frame(akt1BinC.ls), degree=2)

prdMARS <- predict(clsMARS,akt1BinC.ts,type="class",thresh=0.5)

t2=table(prdMARS,Trait.ts)
#       Trait.ts
#prdMARS FALSE TRUE
#  FALSE    54   41
#  TRUE      2    3
t2[2,2]/sum(t2[,2])
t2[1,1]/sum(t2[,1])

prdMARS <- predict(clsMARS,akt1BinC.ts,type="class",thresh=.442)
table(prdMARS,Trait.ts)
prdMARS <- predict(clsMARS,akt1BinC.ts,type="class",thresh=.443)
table(prdMARS,Trait.ts)

# so very sensitive to the threshold, maybe try different degree?

clsMARS <- earth(factor(Trait.ls)~., data=data.frame(akt1BinC.ls), degree=1)

# same problem

# to get ROC curve we will vary the threshold

M=200
specMARS=rep(NA,M)
sensMARS=rep(NA,M)
threshVec=rep(NA,M)
thr=seq(0,1,len=M)
for(i in 1:M){
  prdMARS <- predict(clsMARS,akt1BinC.ts,type="class",thresh=thr[i])
  t2=table(prdMARS,Trait.ts)
  sensMARS[i] <- t2[2,2]/sum(t2[,2])
  specMARS[i] <- t2[1,1]/sum(t2[,1])
}

plot(specMARS,sensMARS,type="l")

#hmm

specMARS
thr[78]
thr[144]

M=200
specMARS=rep(NA,M)
sensMARS=rep(NA,M)
threshVec=rep(NA,M)
thr=seq(.38,.72,len=M)
for(i in 1:M){
  prdMARS <- predict(clsMARS,akt1BinC.ts,type="class",thresh=thr[i])
  t2=table(prdMARS,Trait.ts)
  sensMARS[i] <- t2[2,2]/sum(t2[,2])
  specMARS[i] <- t2[1,1]/sum(t2[,1])
}

specMARS

> thr[33:34]
[1] 0.4346734 0.4363819

############################
## support vector machines

install.packages("kernlab")
library("kernlab")

clsSVM <- ksvm(factor(Trait.ls)~., data=data.frame(akt1BinC.ls), kernel="rbfdot", kpar="automatic", C=60, cross=3, prob.model=TRUE)

prdSVM <- predict(clsSVM,akt1BinC.ts)

t3=table(prdSVM,Trait.ts)
#       Trait.ts
#prdSVM  FALSE TRUE
#  FALSE    38   31
#  TRUE     18   13
t3[2,2]/sum(t3[,2])
t3[1,1]/sum(t3[,1])

# but most would recommend "tuning" the parameters, i.e. change at least C and examine the CV error estimate

clsSVM <- ksvm(factor(Trait.ls)~., data=data.frame(akt1BinC.ls), kernel="rbfdot", kpar="automatic", C=10, cross=3, prob.model=TRUE)

# for instance this does well by that measure

clsSVM <- ksvm(factor(Trait.ls)~., data=data.frame(akt1BinC.ls), kernel="rbfdot", kpar="automatic", C=.1, cross=3, prob.model=TRUE)

> t3[2,2]/sum(t3[,2])
#[1] 0.2954545
> t3[1,1]/sum(t3[,1])
#[1] 0.6785714

# but gives exactly the same out of sample performance

prdSVMp <- predict(clsSVM,akt1BinC.ts,type="probabilities")

svmROC <- roc(Trait.ts,prdSVMp[,1],plot=TRUE)

can compare the methods in terms of confidence intervals for the AUC

ci(svmROC)
ci(rfROC)

# can also test for differences

roc.test(rfROC,svmROC,paired=TRUE)
