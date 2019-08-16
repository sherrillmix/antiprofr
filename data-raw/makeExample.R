antibodyData<-read.csv('example.csv',stringsAsFactors=FALSE,header=FALSE)
antibodyData<-antibodyData[,!apply(is.na(antibodyData),2,all)]
save(antibodyData,file='../data/antibodyData.RData')
tools::resaveRdaFiles('../data/antibodyData.RData')
