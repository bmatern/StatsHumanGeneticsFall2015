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
