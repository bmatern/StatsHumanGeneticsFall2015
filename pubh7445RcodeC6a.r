# Example 6.2
# Reading in Virco data:
virco <- 
read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", 
header=T, sep=",")

# multiple regression

attach(virco)
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
Trait <- as.factor(IDV.Fold > NFV.Fold)

dim(VircoGeno)

# univariate associations

res1 <- rep(NA,99)
for(i in 1:99) res1[i] <- summary(lm(as.numeric(Trait)~as.numeric(VircoGeno[,i])))$coef[2,4]

# gives an error

table(VircoGeno[,i])

# better way

res1 <- rep(NA,99)
for(i in 1:99){
  m1 <- lm(as.numeric(Trait)~as.numeric(VircoGeno[,i]))
  if(is.na(m1$coef[2])==FALSE) res1[i] <- summary(m1)$coef[2,4]
}

summary(res1)
hist(res1)

which(is.na(res1)==TRUE)
which(res1==NA)
which(res1=="NA")

# multiple regression models: problems

m1 <- lm(as.numeric(Trait)~data.matrix(VircoGeno[,1:2]))

table(VircoGeno[,1],VircoGeno[,2])

assocTab <- matrix(NA,99,99)
for(i in 1:99){
  for(j in c(i+1):99){
    assocTab[i,j] <- fisher.test(table(VircoGeno[,i],VircoGeno[,j]))$est 
  }
}

idx1=which(is.na(res1)==FALSE)

assocTab <- matrix(NA,99,99)
for(i in idx1){
  for(j in idx1){
    assocTab[i,j] <- fisher.test(table(VircoGeno[,i],VircoGeno[,j]))$est 
  }
}

assocTab[1,]
table(VircoGeno[,1],VircoGeno[,86])
table(VircoGeno[,1],VircoGeno[,85])
(739*33)/(58*236)

idx2 <- c(1,which(!is.na(assocTab[1,]) & assocTab[1,]<Inf))
m1 <- lm(as.numeric(Trait)~data.matrix(VircoGeno[,idx2]))

sort(summary(m1)$coef[,4])

fisher.test(table(Trait,VircoGeno[,76]))

# compare to previous histogram

hist(summary(m1)$coef[,4])

# making a prediction

sum(m1$coef*c(1,as.numeric(VircoGeno[1,idx2])))

# same as

predict(m1)[1]

plot(as.numeric(Trait[!is.na(Trait)])-1,predict(m1))

# logistic regression

summary(as.numeric(Trait))
m2 <- glm(c(as.numeric(Trait)-1)~data.matrix(VircoGeno[,idx2]),family=binomial)

# smaller data set

m2 <- lm(as.numeric(Trait[1:200])~data.matrix(VircoGeno[1:200,idx2]))

# even smaller: problems

m2 <- lm(as.numeric(Trait[1:50])~data.matrix(VircoGeno[1:50,idx2]))

# Example 6.2 (Creating a classification tree):
attach(virco)
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")
Trait <- as.factor(IDV.Fold > NFV.Fold)
library(rpart)
ClassTree <- rpart(Trait~., method="class", data=VircoGeno)
ClassTree
plot(ClassTree)
text(ClassTree)
rpart(Trait~., method="class", parms=list(split='information'), 
data=VircoGeno)
rpart(Trait~., method="class", parms=list(split='gini'), 
control=rpart.control(minsplit=150, minbucket=50), data=VircoGeno)

# Example 6.3
# Reading in Virco data:
virco <- 
read.csv("http://people.umass.edu/foulkes/asg/data/Virco_data.csv", 
header=T, sep=",")

# Necessary code from Example 6.2:
attach(virco)
VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"]!="-")

# Example 6.3 (Generating a regression tree):
library(rpart)
Trait <- NFV.Fold - IDV.Fold
Tree <- rpart(Trait~., method="anova", data=VircoGeno)
Tree

# Example 6.4
# Reading in FAMuSS data:
fms <- read.delim("http://people.umass.edu/foulkes/asg/data/FMS_data.txt", 
header=T, sep="\t")

# Example 6.4 (Categorical and ordinal predictors in a tree):
attach(fms)
Trait <- NDRM.CH
library(rpart)
RegTree <- rpart(Trait~resistin_c30t+resistin_c398t+
	resistin_g540a+resistin_c980g+resistin_c180g+
	resistin_a537c, method="anova")
RegTree
RegTreeOr <- rpart(Trait~as.numeric(resistin_c30t)+
	as.numeric(resistin_c398t)+as.numeric(resistin_g540a)+
	as.numeric(resistin_c980g)+as.numeric(resistin_c180g)+
	as.numeric(resistin_a537c), method="anova")
RegTreeOr

VircoGeno <- data.frame(virco[,substr(names(virco),1,1)=="P"])
Tree <- rpart(APV.Fold~.,data=VircoGeno)
Tree
plotcp(Tree)
printcp(Tree)
pruneTree <- prune(Tree, cp=0.03)
pruneTree <- prune(Tree, cp=0.07)


