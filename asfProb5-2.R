library(haplo.stats)
hgdpURL <- "http://people.umass.edu/foulkes/asg/data/HGDP_AKT1.txt"
hgdp <- read.delim(file=hgdpURL)

attach(hgdp)


geno=cbind(substr(AKT1.C6024T,1,1), substr(AKT1.C6024T,2,2),
  substr(AKT1.C0756A,1,1), substr(AKT1.C0756A,2,2),
  substr(AKT1.G2375A,1,1), substr(AKT1.G2375A,2,2),
  substr(AKT1.G2347T,1,1), substr(AKT1.G2347T,2,2))
SNPnames <- c("C6024T","C0756A","G2375A", "G2347T")

table(Population)
dim(geno)

#This is for the first observation.  We'll do the same thing in a loop.
#Im not sure what's going on here lol.
i=1
id=levels(Population)[i]
dat=geno[Population==id,]
h1=haplo.em(dat, locus.label=SNPnames,
  control=haplo.em.control(min.posterior=1.0e-4))
h1
names(h1)
h1$hap.prob
h1$haplotype
#maxNoHap is an upper bound for matrix size.
maxNoHap=40
hapProb=matrix(NA, 52, maxNoHap)
hapProb
#this is a vector of lists?  I think.
hapMat=vector("list",52)
for(i in 1:52){
  id=levels(Population)[i]
  dat=geno[Population==id,]
  h1=haplo.em(dat, locus.label=SNPnames,
    control=haplo.em.control(min.posterior=1.0e-4))
  nh=dim(h1$haplotype)[1]
  hapProb[i,1:nh]=h1$hap.prob
  hapMat[[i]]=h1$haplotype
}

dat=numeric(0)
dat

for(i in 1:52){
  #nh = number of haplotypes
  nh=dim(hapMat[[i]])[1]   
  for(j in 1:nh) dat=c(dat,hapMat[[i]][j,])
}

#now dat is a 620 vector of nucleotides for haplotypes.
dat
#This is the number of haplotypes found.  each haplotype is a series of 4 nucleotides.  We're gonna paste them togesther.
#PASTe seems to be the same as concat.  
length(dat)/4
#converting the vector into a matrix
ddat1=matrix(dat, byrow=T, ncol=4)
ddat1
dat2=rep(NA, 155)
for(i in 1:155) dat2[i]=paste(ddat1[i,1],ddat1[i,2],ddat1[i,3],ddat1[i,4],sep="")
table(dat2)
hapProb[1:5,1:10]
#t = transpose.
t(hapProb[1:5,1:10])
#c seems to turn a matrix into a vector of some sort.
#c function puls vallues vertically first, that's why he transposed it first.
c(t(hapProb[1:5,1:10]))
hapProbPop=c(t(hapProb))[!is.na(c(t(hapProb)))]
length(hapProbPop)

#We're doing the weird example then loop thing again?
#U tgubj were trying to make a vector of histograms
i=1
names(table(dat2))[i]
hid=names(table(dat2))[i]
hapProbPop[dat2==hid]
#par controls grapical parameters.  mfrow is a 3x3 series of histograms.  Ok that's neat.
par(mfrow=c(3,3))
for(i in 1:9){
  hid=names(table(dat2))[i]
  hist(hapProbPop[dat2==hid])
}
#mar = margins.  Adjusting space between the histograms.
#
par(mar=c(2,2,2,1))
for(i in 1:9){
  hid=names(table(dat2))[i]
  hist(hapProbPop[dat2==hid], xlab="", main=paste("hap:", hid))
}



detach(hgdp)