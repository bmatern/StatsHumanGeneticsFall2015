######################################################################
# R commands      Statistics for Human Genetics and Molecular Biology
# Lab 2 2015                                    University of Minnesota
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


##############################
# R as a calculator
##############################

2+3
2^3
2+3-2^3
(2+3-2^3)*4
(2+3-2^3)/4
2+3^4
1:4
2+(1:4)^4
log(1:4)
log2(1:4)
log(1:4, base=3)
sin(0.5)
7%% 3  ## answer 1, modulo reduction
abs(-3)
exp(1)
sqrt(3)
x<-3
x
x<-4
x

##############################
# Creating simple vectors
##############################

x <- c(1,3.5,-28.4,10) #numerical vector
x
y<-c("cat","dog","mouse","monkey") #character
z<-c(TRUE,TRUE,TRUE,FALSE,FALSE) #logical vector
x<-1:10
seq(1,10)
seq(3, 9, by=3)
rep(2, 10)
rep(c(1,2,3),5)
log(seq(1,2, by=0.1))
x<-c(1, 5, 10, NA, 15)
sum(x)
sum(x, na.rm=T)
prod(x, na.rm=T)
mean(x, na.rm=T)
x<-1:10
cumsum(x) ## running sum
cummax(x) ## running maximum
cummin(x) ## running minimum

###################################
# Accessing Elements in a Vector
###################################
y <- c(18,32,15,-7,12,19)
length(y)
y[3:5] ##position in vector as positive integer
y[-c(1,5,6)] ## exclude: use negative integers
y<15 
y[y<15]
which(y==32)
x<-1:10
match(y,x)
colors<-c("red", "blue", "pink")
which(colors=="yellow")
x<-c(1,5,10, NA, 15)
which(is.na(x))
which(!is.na(x))

##############################
# Factors
# vector with categories
##############################
colors<-c(1,1,2,3)
colors<-factor(colors, label=c("red", "green", "blue"))
table(colors)



##############################
# Matrices
##############################
help(cbind)
y<-c(8,32, 15, -7, 2, 19)
x<-1:6
mat<-cbind(x,y)
help(rbind)
dim(mat)  ## check dimension
ncol(mat) ## the number of columns of a matrix
nrow(mat) ## the number of rows of a matrix
mat[2,3] # the value in the 2nd row and the 3rd column
mat[1:3,]  ## the first three row of mat
mat[,2]  ## the 2nd column of mat
mat[-1,] ## exclude the first row
newmat<-matrix(1:9, nrow=3) ## create new matrix
newmat
rowMeans(newmat)
colMeans(newmat)
m<-matrix(1:9, nrow=3, byrow=T) ## fill row first
colnames(m)<-c("a", "b", "c")
rownames(m)<-c("r1", "r2", "r3")
vect<-as.vector(newmat)

##############################
# Matrices Multiplication
##############################
mat<-matrix(1:9, nrow=3)
mat^2
mat%*%mat

##############################
# Arrays
##############################

myarray<-array(1:64, dim=c(4,4,4))
myarray
myarray[1,2,3]



##############################
# Data Frames
##############################
muscle <- rnorm(n=10,mean=3,sd=1)
sex <- factor(rep(c("M","F"),c(6,4)))
speed  <- rep(0,10)
speed[1:6] <- rnorm(6,30-2*muscle[1:6],2)
speed[7:10] <- rnorm(4,40-2*muscle[7:10],2)
mydata <- data.frame(y=speed,x1=muscle,x2=sex)
summary(mydata)
str(mydata)
temp <- lm(y~x1+x2,data=mydata)  
summary(temp)

##############################
# Lists
##############################
x <- list(one=c(18:36),two=c("AK","AL","AZ"),three=c(T,T,F,T),four=matrix(1:12,3,4))
x
x[[1]][3:6]
x$one[3:6]
y<-unlist(x)
str(y)


##############################
# For loops
##############################

for (i in 1:100){
	d<-Sys.time()
	print(paste("Now is", d, sep=" "))
	print(i*i)
}


##############################
# Generate random number 
##############################

## from uniform distribution 

x<-runif(100)
hist(x)

## from normal distribution 

y<-rnorm(100)
plot(density(y))



##################
# End of lab1.R
##################
