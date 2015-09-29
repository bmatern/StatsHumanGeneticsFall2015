######################################################################
# R commands      Statistics for Human Genetics and Molecular Biology
# Lab 3  2015                                 University of Minnesota
######################################################################
# This file contains the R commands for the lab.
#
# Lines beginning with the symbol '#' are comments in R.  All other
# lines contain code.
#
# In R for Windows, you may wish to open this file from the menu bar
# (File:Display file); you can then copy commands into the command
# window.  (Use the mouse to highlight one or more lines; then
# right-click and select "Paste to console".)
######################################################################


##########################################################################
# Import data 
#############################################################################
## Need to creat a dataset with missing values in excel and save it as test.csv (comma delimited) for this exercise
getwd()
data<-read.csv("test.csv", header=T)
str(data)
rm<-is.na(data$Values)
rm<-which(is.na(data$Values))
newdata<-data[-rm,]
######Alternative using !is.na
keep<-which(!is.na(data$Values))
newdata<-data[keep,]
##############################
## Export data 
##############################
write.table(newdata, file="test2.txt", row.names=F, quote=F, sep="\t")

############# FAMuss example 
getwd() ## change it using setwd()
fmsURL<-"http://people.umass.edu/foulkes/asg/data/FMS_data.txt"
fms<-read.delim(file=fmsURL, header=TRUE, sep="\t")
colnames(fms)
dim(fms)  ## check the dimension of the data
str(fms[,1:10]) ## check the structure of the data 
fms$id[1:10] 
fms[1,1:10]
fms$pre.BMI
fms$actn3_rs540874
index<-match(c("actn3_rs540874","pre.BMI"), colnames(fms)) 
dat<-fms[!is.na(fms$pre.BMI) & !is.na(fms$actn3_rs540874) , index] ## observations without NA
attach(fms)
mean(pre.BMI, na.rm=T)
dURL<-"http://itl.nist.gov/div898/strd/univ/data/PiDigits.dat"
pidigits<-read.table(file=dURL, skip=60)
table(pidigits)
str(pidigits)
prop<-table(pidigits)/5000 ##proportions
prop
barplot(prop, xlab="digit", ylab="proportion")
abline(h=0.1, lty=2)
##scan
x<-scan()
str(x)

##############################
## &, |
##############################
x<-c(T, T, F, F)
y<-c(T, F, T, F)
mat<-cbind(x,y)
mat
and<- x & y
or<-x | y
and 
or


##############################
# Export data 
##############################
newdata<-fms[,1:3]
write.table(newdata, file="newdata.txt", sep="\t")
save(newdata, file="newdata.RData")


##############################
# Exercise 1
##############################
# 1. create a file using excel 
# 2. save it as csv
# 3. read it into R using read.csv
##############################
# Exercise 2
##############################
# 1. create a file using excel 
# 2. save it as txt (tab-delimited)
# 3. read it into R using read.table
##############################


##############################
# Generate Random Number  
###############################
##################################################
# generate 100 points from Normal(mu=0, sigma=2)
##################################################
x<-rnorm(100, mean=0, sd=2)
mean(x)
hist(x)
plot(density(x))
mean(x)
quantile(x)
##################################################
# generate 100 points from uniform distribution(0,1)
##################################################
y<-runif(100)
sum(y > 0.5)

##################################################
# Toss a fair coin five times 
##################################################

h<-rbinom(1, 5, p=0.5)

##################################################
# sample from a vector
##################################################
x<-1:10
samp<-sample(x)
samp2<-sample(x, 5)
samp
samp2

##############################
# Graphics 
##############################
library(MASS)
head(mammals)
str(mammals)
br<-seq(0, max(mammals$body)+100, by=100)
hist(mammals$body, freq=T, breaks=br, xlab="Body Size (lbs)", main="Distribution body size of mammals")
m<-median(mammals$body)
size<-ifelse(mammals$body > m, "large", "small")
mammals$size<-size
head(mammals)
plot(mammals$body, mammals$brain, pch=16, xlab="Body Size (lbs)", ylab="Brain Size (g)")
rm<-which(mammals$body > 1000)
mamm<-mammals[-rm,]
dim(mamm)
plot(mamm$body, mamm$brain, pch=16, xlab="Body Size (lbs)", ylab="Brain Size (g)")
small<-mammals[mammals$size=="small",]
plot(small$body, small$brain, pch=16)

############################
x<-rnorm(100)
y<-3 + 2*x + rnorm(100)
plot(x, y, type="p", pch=16)
## plot normal(0,1) distribution
plot(x<-seq(-4,4, by=0.1), dnorm(x, mean=0, sd=1), type="l", xlab="X", ylab="Density", main="Normal distribution")
hist(x)
pdf("myplot.pdf")
plot(x<-seq(-4,4, by=0.1), dnorm(x, mean=0, sd=1), type="l", xlab="X", ylab="Density", main="Normal distribution")
dev.off()

####################
# Installing packages
#  
######################
source("http://www.bioconductor.org/biocLite.R")
biocLite("UsingR")
library("UsingR")


### generate 1000 points from chi square distribution
x<-rchisq(1000, df=1)
quantile(x)
maint<-expression(paste(chi^2, " distribution", sep=""))
plot(density(x), main=maint, xlab="x")
simple.hist.and.boxplot(x)



##############################
# for loops & ifelse
##############################

for (i in 1:10){
	if (i > 5){
		print(i)
	}
}
## if .. else 
x<-rnorm(10, mean=5, sd=1)

for (i in x){
	if (i < 5){
		msg<-paste(i, " is equal or less than 5", sep="")
		print(msg)
	} else {
		msg<-paste(i, " is larger than 5", sep="")
		print(msg)
	}
	
}


str(fms[,1:10])
attach(fms)
table(actn3_rs540874)
geno<-ifelse(actn3_rs540874=="GG", 1, 0)
table(geno)
geno2<-as.numeric(actn3_rs540874)
table(geno2, useNA="always")
x<-rnorm(100)
xx<-ifelse(x > 0, 1, 0)
table(xx)



##############################
# Functions
##############################
x<-rnorm(100)
mysd<-function(x){
	ans<-sqrt(var(x))
	return(ans)
}
mysd(x)
sd(x)

myfun<-function(x){
	ans<-x+3
	return(ans)
}
myfun(5)

##############################
# Likelihood
##############################

likelihood<-function(p, n=40, k=25){
	ans<-dbinom(x=k, size=n, prob=p)
	return(ans)
}

getwd()
pdf("MLE.pdf")
plot(x<-seq(0, 1, by=0.01), likelihood(x, n=40, k=25), type="l", xlab="p", ylab="Likelihood", main="Toss a coin 40 times and get 25 heads")
abline(v=0.625,col="grey", lty=2)
text(0.3, 0.12, "Maximum Likelihood p=0.625", bty="n")
text(0.3, 0.115, "(the most likely value of p)", bty="n")
dev.off()


##############################
# apply, sapply
##############################
m<-matrix(1:20, nrow=10, ncol=2) ## create a matrix of 10 rows x 2 columns
apply(m, 1, mean) ## mean of the rows: rowMeans
apply(m,2, mean)  ## mean of the columns: colMeans 
x<-1:10
sapply(x, myfun)
x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
sapply(x, quantile)

###########################
# Exercise: Tossing a coin a lot of times 
###########################



##################
# End of lab2.R
##################




