# 1)
# Foulkes 6.2
virco <-  read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv",   header=T, sep=",")

attach(virco)
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
VircoGeno
Trait <- as.factor(APV.Fold > IDV.Fold)
library(rpart)
ClassTree <- rpart(Trait~., method="class", data=VircoGeno)
ClassTree
plot(ClassTree)
text(ClassTree)

pruneTree <- prune(ClassTree,cp=.03)

pruneTree
plot(pruneTree)
text(pruneTree)

rpart(Trait~., method="class", parms=list(split='information'), 
      data=VircoGeno)
rpart(Trait~., method="class", parms=list(split='gini'), 
      control=rpart.control(minsplit=150, minbucket=50), data=VircoGeno)
detach(virco)

# 2)
# Foulkes 6.5
#fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt",  header=T, sep="\t")
#virco <-  read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv",   header=T, sep=",")
#attach(virco)
#attach(fms)
#fms
#VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
#SNPnames <- c("resistin_a537c", "resistin_c180g", "resistin_c30t"
#              , "resistin_c398t" , "resistin_c980g",   "resistin_g540a"
#              , "actn3_1671064" , "actn3_r577x"         
#              , "actn3_rs1815739" , "actn3_rs540874"  )
#sort(colnames(fms))
#FmsGeno <- data.frame(fms[,SNPnames])

#library(rpart)
#Trait <- NFV.Fold - IDV.Fold
#Trait <- NDRM.CH
#Tree <- rpart(Trait~., method="anova", data=VircoGeno)
#Tree <- rpart(Trait~., method="anova", data=FmsGeno)
#Tree
#detach(virco)

fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt",  header=T, sep="\t")
attach(fms)
Trait <- NDRM.CH
library(rpart)
RegTree <- rpart(Trait~resistin_c30t+resistin_c398t+
                   resistin_g540a+resistin_c980g+resistin_c180g+
                   resistin_a537c+actn3_1671064+actn3_r577x+actn3_rs1815739+
                    actn3_rs540874 , method="anova")
RegTree
RegTreeOr <- rpart(Trait~as.numeric(resistin_c30t)+
                     as.numeric(resistin_c398t)+as.numeric(resistin_g540a)+
                     as.numeric(resistin_c980g)+as.numeric(resistin_c180g)+
                     as.numeric(resistin_a537c)+
                     as.numeric(resistin_a537c)+
                     as.numeric(actn3_1671064)+
                     as.numeric(actn3_r577x)+
                     as.numeric(actn3_rs1815739)+
                     as.numeric(actn3_rs540874)
                   
                   , method="anova")
RegTreeOr

#VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"])
#Tree <- rpart(APV.Fold~.,data=VircoGeno)
#Tree
#plotcp(Tree)
#printcp(Tree)
#pruneTree <- prune(Tree, cp=0.03)
#pruneTree <- prune(Tree, cp=0.07)
















# 3)
# Foulkes 7.1
# Fit a random forest to the Virco data to determine 
# which protease mutations are most highly associated 
# with SQV fold resistance as measured by SQV.fold.
virco <-  read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", header=T, sep=",")

#install.packages("randomForest")
library(randomForest)
attach(virco)
#Trait <- NFV.Fold - IDV.Fold
Trait <- SQV.Fold
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
Trait.c <- Trait[!is.na(Trait)]
VircoGeno.c <- VircoGeno[!is.na(Trait),]
RegRF <- randomForest(VircoGeno.c, Trait.c, importance=TRUE)
RegRF
varImpPlot(RegRF,main="")


# 4)
# Foulkes 7.2
# Apply logic regression to the Virco data to characterize, 
# with a logic structure, the association between protease mutations 
# and SQV fold resistance as measured by SQV.fold.
virco <-  read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv",     header=T, sep=",")
attach(virco)
Trait <- SQV.Fold
#Trait <- NFV.Fold - IDV.Fold
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
Trait.c <- Trait[!is.na(Trait)]
VircoGeno.c <- VircoGeno[!is.na(Trait),]
#install.packages("LogicReg")
library(LogicReg)
VircoLogicReg <- logreg(resp=Trait.c, bin=VircoGeno.c, select=1)
plot(VircoLogicReg)
VircoLogicReg
VircoLogicRegMult <- logreg(resp=Trait.c, bin=VircoGeno.c, select=2, 
                            ntrees=2, nleaves=8)
plot(VircoLogicRegMult)
VircoLogicRegMult
detach(virco)

# 5) 
# Foulkes 7.4 
# 7.4. Apply a multivariable adaptive regression spline to address 
# the question of Problem 7.3. 
# Compare and contrast your findings.

# 7.3
#Apply Monte Carlo logic regression to the FAMuSS data to characterize 
#the importance of SNPs in the actn3 and resistin genes in predicting 
#change in non-dominant arm muscle strength as measured by NDRM.CH.

#virco <- read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv",     header=T, sep=",")
# Example 7.6 (Monte Carolo logic regression):
library(LogicReg)
attach(virco)

fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt",  header=T, sep="\t")
attach(fms)

#Trait <- SQV.Fold
Trait <- NDRM.CH

#VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
FmsGeno <- data.frame(!is.na(fms[,c("resistin_a537c", "resistin_c180g", "resistin_c30t"
                             , "resistin_c398t" , "resistin_c980g",   "resistin_g540a"
                              , "actn3_1671064" , "actn3_r577x"         
                               , "actn3_rs1815739" , "actn3_rs540874" 
                             )]))

#FmsGeno <- is.na(FmsGeno)

Trait.c <- Trait[!is.na(Trait)]
Trait
Trait.c
#VircoGeno.c <- VircoGeno[!is.na(Trait),]
FmsGeno.c <- FmsGeno[!is.na(Trait),]
#VircoLogicRegMCMC <- logreg(resp=Trait.c, bin=VircoGeno.c, select=7)
FmsLogicRegMCMC <- logreg(resp=Trait.c, bin=FmsGeno.c, select=7)
#plot(sort(VircoLogicRegMCMC$single), xlab="Sorted SNPs", ylab="Number of 
#     selected models")
plot(sort(FmsLogicRegMCMC$single), xlab="Sorted SNPs", ylab="Number of selected models")
#names(VircoGeno)[order(VircoLogicRegMCMC$single)]
names(FmsGeno)[order(FmsLogicRegMCMC$single)]


#virco <- read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv",   header=T, sep=",")


#attach(virco)
#Trait <- NFV.Fold - IDV.Fold
#VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
#Trait.c <- Trait[!is.na(Trait)]
#VircoGeno.c <- VircoGeno[!is.na(Trait),]
#install.packages("earth")
library(earth)
#VircoMARS <- earth(Trait.c~., data=VircoGeno.c, degree=2)
FmsMARS <- earth(Trait.c~., data=FmsGeno.c, degree=32)
FmsMARS
summary(FmsMARS)
evimp(FmsMARS)

detach(virco)
