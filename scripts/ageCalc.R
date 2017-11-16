#Load in the libraries we will need
library(metafor)

#Move into script directory
setwd(dirname(sys.frame(1)$ofile))

#Load in  meta analysis data
metData = read.csv('../data/metaAvData.csv')

#Get age vector. Remove missing data and studies that only reported a range
metAge = metData$Age[is.na(metData$Age)==FALSE]
metAge = as.numeric(as.character(metAge[-grep('-',metAge)]))

#Show mean age to user
print(paste('Mean Age:',round(mean(metAge),2)))
print(paste('SD:',round(sd(metAge),2)))