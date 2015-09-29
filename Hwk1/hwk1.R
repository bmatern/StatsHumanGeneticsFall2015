#Homework 1
#Ben Matern

#1)
#2.1
miles = c(65311, 65624, 65908, 66219, 66499, 66821, 67145, 67447)
x = diff(miles)
max(x)
mean(x)
min(x)
#2.2
commute = c(17, 16, 20, 24, 22, 15, 21, 15, 17, 22)
max(commute)
mean(commute)
min(commute)
commute = replace(commute, commute==24, 18)
commute
sum(commute >= 20)
sum(commute < 17) / length(commute)
#2.3
bill = c(46, 33, 39, 37, 46, 30, 48, 32, 49, 35, 30, 48)
sum(bill)
min(bill)
length(bill[bill>40])
length(bill[bill>40])/length(bill)

#2)
#2.5
x = c(1,3,5,7,9)
y = c(2,3,5,7,11,13)
x+y
y[x]
y[y>=7]

#3)
#2.6
x = c(1,8,2,6,3,8,5,5,5,5)
sum(x)/10
log10(x)
log(x, 10)
(x-4.4)/2.875
max(x) - m

#4)
c(rep(1,5), rep(2,4), rep(3,3), rep(4,2), 5)

#5)
mat <- matrix(c(1:13,NA,14:15),nrow=4)
mat
badIndices = !is.na(mat[,4]) & mat[,4] < 14
badIndices  
result = mat[-badIndices,]
result

#6)
fmsURL = "http://people.umass.edu/foulkes/asg/data/FMS_data.txt"
fms <- read.delim(file=fmsURL, header=T, sep="\t")
fms
#a
gadgets <- 1:600
thingies <- 301:900
objects  <- 401:1000
superobjects <- objects * 3
#b
myNewData <- data.frame(y=gadgets,x1=thingies,x2=objects,x3=superobjects)
str(myNewData)
write.table(myNewData, file="/Users/bmatern/school/Fall2015/Stats4HumGenomics/Hwk1/benoutput.txt", row.names=F, quote=F, sep="\t")
reloadedData <- read.csv("/Users/bmatern/school/Fall2015/Stats4HumGenomics/Hwk1/benoutput.txt",sep="\t")
str(reloadedData)
#c
save(myNewData, file="/Users/bmatern/school/Fall2015/Stats4HumGenomics/Hwk1/newdata.RData")
load("/Users/bmatern/school/Fall2015/Stats4HumGenomics/Hwk1/newdata.RData")
str(myNewData)

#8
#a
#10,000,000x10 matrix
mat <- matrix(runif(100000000, min=0, max=1),ncol=10)
beginTime <- proc.time()
for(i in 1:10)
{
  meanValue <- mean(mat[,i])
}
endTimeMean <- proc.time() - beginTime



#b
#10,000,000x10 matrix
mat <- matrix(runif(100000000, min=0, max=1),ncol=10)
beginTime <- proc.time()
for(i in 1:10)
{
  means <- apply(mat,2,mean)
}
endTimeApply <- proc.time() - beginTime

#9
#a
nucs<-c("A","G","C","T")
dna<-sample(nucs, 10000, replace=T)
#b
table(dna)
#c

contigTable<-function(nucs, dna)
{
  results <- matrix(rep("",128),nrow=2)
  resultsIndex = 1
  for(i in 1:4)
  {
    charI <- nucs[i]
    for(j in 1:4)
    {
      charJ <- nucs[j]
      for(k in 1:4)
      {
        charK <- nucs[k]
        sequence <- paste(charI, charJ, charK, collapse="")
        occuranceCount <- 0
        
        for(x in 1:(10000-3))
        {
          if(dna[x]==charI && dna[x+1]==charJ && dna[x+2]==charK)
          {
            occuranceCount <- 1 + occuranceCount
          }
        }
        results[1,resultsIndex] <- sequence
        results[2,resultsIndex] <- occuranceCount
        resultsIndex <- resultsIndex + 1
      }
    }
  }
  
  return(results)
}
myTable = contigTable(nucs,dna)
myTable

#10
factorialFunction<-function(upperValue)
{
  if(upperValue < 1)
  {
    return(0)
  }
  results <- 1
  for(x in 2:upperValue)
  {
    results <- results * x
  }
  return(results)
}

factorialFunction(10)


