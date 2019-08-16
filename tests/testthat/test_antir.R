
context("taxa functions")
test_that("Test read.anti",{
  tmp<-tempfile()
  #check if data.frame and file read in the same
  write.table(antibodyData,tmp,row.names=FALSE,col.names=FALSE,sep=',')
  expect_equal(readAnti(antibodyData),readAnti(tmp))
  expect_message(readAnti(antibodyData),'32 plate.*8 row')
  expect_message(readAnti(antibodyData,vocal=FALSE),NA)
  tmp<-antibodyData
  tmp[,1]<-dnar::fillDown(tmp[,1])
  expect_equal(readAnti(antibodyData),readAnti(tmp))
  tmp2<-tmp
  tmp2[,]<-NA
  tmp<-rbind(antibodyData,tmp2)[order(rep(1:nrow(antibodyData),2)),]
  expect_equal(readAnti(antibodyData),readAnti(tmp))
  tmp2[,2]<-''
  tmp<-rbind(antibodyData,tmp2)[order(rep(1:nrow(antibodyData),2)),]
  expect_equal(readAnti(antibodyData),readAnti(tmp))
  tmp<-cbind(tmp,NA,NA,NA,NA)
  expect_equal(readAnti(antibodyData),readAnti(tmp))
  expect_error(readAnti(antibodyData[1:14,]),'rows')
  expect_error(readAnti(antibodyData[1:16,]),NA)
  expect_error(readAnti(antibodyData[1:14,],nrow=7),NA)
})
