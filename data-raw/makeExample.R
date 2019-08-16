antibodyData<-read.csv('example.csv',stringsAsFactors=FALSE)
save(antibodyData,file='../data/antibodyData.RData')
tools::resaveRdaFiles('../data/antibodyData.RData')
