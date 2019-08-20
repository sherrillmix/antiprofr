#' Read antibody profiling data
#'
#' Read in antibody profiling data from 96 well plates using positive and negative controls and replicate serial dilutions of patient plasma. Negative control data is assumed to be in the first 4 columns of the first row of data for each plate. The mean of the negative contol is subtracted from all values in the \code{sub} column of the output.. Positive control data with a known antibody is assumed to be in columns 5-8 and 9-12 of the first row.
#'
#' @param dat a file containing the raw data or a data.frame containing raw data. The function assumes that the first column of the first row of each plate contains the patient ID and the remainder of that plate is empty.
#' @param positiveThreshold a single numeric threshold to use for determining a threshold by estimating the OD for this amount of positive control antibody
#' @param dilutions
#' @param p24Dilutions
#' @param vocal
#' @param nrows
#' @return a data.frame with a row for each well in the raw plate data with columns:
#' \itemize{
#'  \item pat The patient identifier from the first column of the raw plate data
#'  \item antigen The antigen identifier from the second column of the raw plate data
#'  \item dilution The serum dilution for the given well
#'  \item pos TRUE if the well was a positive control. Otherwise FALSE
#'  \item neg TRUE if the well was a negative control. Otherwise FALSE
#'  \item od The raw OD value read from the plate data
#'  \item plate The plate number for a well. Plates are numbered starting from 1 for the topmost plate.
#'  \item noAnt The OD value for the negative control well corresponding to this serum dilution
#'  \item sub The raw value with the negative control mean subtracted
#' }
#' @export
#' @seealso \code{\link{calcP24Cut}}, \code{\link{calcCross}}, \code{\link{plotAnti}}
#' @examples
#' fakeDat<-1:8 %*% t(1/3^rep(0:3,3))
#' fakeDat[1,]<-c(rep(0.1,4),50*1/3^rep(0:3,2))
#' raw<-cbind(data.frame('pat'='Test',antigen=LETTERS[1:8]),fakeDat)
#' tmpFile<-tempfile()
#' write.table(raw,tmpFile,row.names=FALSE,col.names=FALSE,sep=',')
#' readAnti(tmpFile)
#' readAnti(raw)
readAnti<-function(dat,positiveThreshold=25,dilutions=300*3^(0:3),p24Dilutions=240/3^(0:3),vocal=TRUE,nrows=8){
  if(!is.data.frame(dat))dat<-utils::read.csv(dat,header=FALSE,stringsAsFactors=FALSE)
  dat<-dat[dat[,2]!=''&!is.na(dat[,2]),]
  if(nrow(dat)%%nrows!=0)stop('Found ',nrow(dat),' rows. Not evenly divided by nrows=',nrows)
  patLines<-seq(1,nrow(dat),8)
  if(vocal)message('Reading in ',length(patLines),' plates of ',nrows,' rows')
  stacked<-do.call(rbind,lapply(1:length(patLines),function(ii){
    xx<-patLines[ii]
    pat<-trimws(dat[xx,1])
    noAnt<-unlist(dat[xx,c(2+1:4)])
    p24Ab<-unlist(dat[xx,c(2+5:12)])
    raws<-unlist(dat[xx+1:7,2+1:12])
    ants<-trimws(dat[xx+1:7,2])
    names(noAnt)<-dilutions
    out<-data.frame('pat'=pat,'antigen'=rep(ants,12),'dilution'=rep(rep(dilutions,3),each=7),'pos'=FALSE,'neg'=FALSE,'od'=raws,stringsAsFactors=FALSE)
    ns<-c(length(noAnt),length(p24Ab))
    out<-rbind(out,data.frame('pat'=pat,'antigen'=rep(c('noAntigen','p24Control'),ns),'dilution'=c(dilutions,rep(p24Dilutions,2)),'pos'=rep(c(F,T),ns),'neg'=rep(c(T,F),ns),'od'=c(noAnt,p24Ab),stringsAsFactors=FALSE))
    out$plate<-ii
    out$noAnt<-noAnt[as.character(out$dilution)]
    out$sub<-out$od-out$noAnt
    return(out)
  }))
  p24<-stacked[stacked$pos,]
  cuts<-calcP24Cut(p24$plate,p24$od,p24$dil,positiveThreshold=positiveThreshold)
  stacked$cut<-cuts[as.character(stacked$plate)]
  xx<-stacked[!stacked$pos&!stacked$neg,]
  tmp<-paste(xx$antigen,xx$pat,xx$plate,sep='___')
  crosses<-by(xx[,c('sub','dilution','cut')],tmp,function(yy)calcCross(yy$sub,yy$dil,yy$cut[1]))
  stacked$cross[!stacked$pos&!stacked$neg]<-crosses[tmp]
  stacked
}
calcP24Cut<-function(plate,od,dil,positiveThreshold=25){
  odDil<-data.frame('od'=od,'dil'=dil,'plate'=plate,'logDil'=log(dil))
  out<-sapply(unique(plate),function(ii){
    thisDat<-odDil[odDil$plate==ii,]
    fit<-stats::lm(I(log(od))~logDil,thisDat)
    od25<-exp(stats::predict(fit,data.frame('logDil'=log(positiveThreshold))))
    return(od25)
  })
  names(out)<-unique(plate)
  out
}
calcCross<-function(od,dil,odCut=25){
  odDil<-data.frame('od'=suppressWarnings(log(od)),'dil'=log(dil))
  odDil<-odDil[!is.na(odDil$od)&odDil$od>log(.001),]#&odDil$od>-3,]
  if(length(unique(odDil$dil))<2)return(0)
  fit<-stats::lm(od~dil,data=odDil)
  if(fit$coef['dil']>0)return(0)
  exp(-(fit$coef['(Intercept)']-log(odCut))/fit$coef['dil'])
}
plotAnti<-function(stacked,p24Cut=25,plotBoth=TRUE,lab='Dilution(OD450=p24 25pg)'){
  ylim<-range(stacked$od)
  p24<-stacked[stacked$pos,]
  p24<-p24[order(p24$plate,p24$dilution),]
  negs<-stacked[stacked$neg,]
  negs<-negs[order(negs$plate,negs$dilution),]
  ylimP24<-range(p24$od,negs$od)
  graphics::par(mfrow=c(ceiling(length(unique(p24$plate))/2),2),mar=c(4.5,4,2.5,.4))
  for(ii in unique(p24$plate)){
    thisDat<-p24[p24$plate==ii,]
    graphics::plot(thisDat$dilution,thisDat$od,xlim=range(p24$dilution),ylim=ylimP24,xlab='',ylab='OD450',las=1,main=paste("p24 Standard Plate",ii),log='xy')
    graphics::abline(h=negs[negs$plate==ii,'od'],col='#0000FF22',lty=2)
    for(ii in unique(thisDat$dilution))graphics::segments(ii,min(thisDat[thisDat$dilution==ii,'od']),ii,max(thisDat[thisDat$dilution==ii,'od']),col='red')
    thisDat$logDil<-log(thisDat$dil)
    fit<-stats::lm(I(log(od))~logDil,thisDat)
    fakeDat<-data.frame('logDil'=seq(1e-2,10,.01))
    pred<-stats::predict(fit,fakeDat)
    od25<-exp(stats::predict(fit,data.frame('logDil'=log(p24Cut))))
    graphics::abline(v=p24Cut,h=od25,lty=2)
    graphics::lines(exp(fakeDat$logDil),exp(pred))
  }
  for(ii in unique(stacked$antigen[!stacked$antigen %in% c('noAntigen','p24Control')])){
    graphics::par(mfrow=c(ceiling(length(unique(stacked$pat))/ifelse(plotBoth,1,2)),2),mar=c(4.5,4,2.5,.4))
    for(jj in unique(stacked$pat[stacked$antigen==ii])){
      thisDat<-stacked[stacked$pat==jj&stacked$antigen==ii,]
      thisPlate<-unique(stacked[stacked$pat==jj&stacked$antigen==ii,'plate'])
      thisDat<-thisDat[order(thisDat$dil,thisDat$plate),]
      neg<-tapply(stacked$od[stacked$neg&stacked$pat==jj&stacked$plate==thisPlate],stacked$dilution[stacked$neg&stacked$pat==jj&stacked$plate==thisPlate],mean)
      thisDat$sub<-thisDat$od-neg[as.character(thisDat$dilution)]
      subMean<-tapply(thisDat$sub,thisDat[,c('dilution')],stats::median)
      above5<-stats::approx(subMean,log10(as.numeric(names(subMean))),thisDat$cut[1])$y
      if(plotBoth){
        graphics::plot(1/thisDat$dil,thisDat$od,main=sprintf('%s %s',ii,jj),xlab='Plasma Dilution',log='x',ylab='OD450',las=1,xaxt='n',ylim=ylim,mgp=c(2.5,.7,0))
        dnar::logAxis(1,axisMin=1e-4)
      }
      graphics::plot(1/thisDat$dil,ifelse(thisDat$sub<.01,.01,thisDat$sub),main=sprintf('%s %s\n %s=%0.0f',ii,jj,lab,thisDat$cross[1]),xlab='Plasma Dilution',log='xy',ylab='OD450 (-negative)',las=1,xaxt='n',ylim=c(.01,ylim[2]),mgp=c(2.5,.7,0),cex=2)
      dnar::logAxis(1,axisMin=1e-4)
      graphics::lines(1/as.numeric(names(subMean)),subMean,lty=3)
      graphics::abline(h=unique(thisDat$cut),lty=2)
      graphics::abline(v=1/unique(thisDat$cross),lty=2,col='red')
      thisDat$logDil<-log(thisDat$dil)
      if(sum(thisDat$sub>.001)>2){
        fit<-stats::lm(I(log(sub))~logDil,thisDat[thisDat$sub>0.001,])
        fakeDat<-data.frame('logDil'=seq(1e-2,10,.01))
        pred<-stats::predict(fit,fakeDat)
        graphics::lines(1/exp(fakeDat$logDil),exp(pred),lty=1)
      }
    }
  }
}
convertStackedToMat<-function(stacked){
  uniq<-stacked[!duplicated(paste(stacked$pat,stacked$antigen,stacked$plate))&!stacked$pos&!stacked$neg,c('pat','antigen','plate','cross')]
  out<-tapply(uniq$cros,list(uniq$pat,uniq$antigen),mean)
  return(out)
}
plotHeat<-function(odMat,filterLess=min(odMat,na.rm=TRUE),scaleMain='Dilution reaching OD450(p24=25pg)',groupings=NULL,...,groupSpace=.5,yLabs=data.frame(colnames(odMat)),filterMore=NULL){
  if(any(selector<-!rownames(odMat)%in%unlist(groupings)))groupings<-c(groupings,list(rownames(odMat)[selector]))
  if(any(table(unlist(groupings)))>1)stop('Rowname appears in more than 1 group')
  group<-rep(1:length(groupings),sapply(groupings,length))
  newGroup<-c(FALSE,group[-1]!=group[-length(group)])
  if(any(!unlist(groupings) %in% rownames(odMat)))stop('Extra group names found')
  odMat<-odMat[unlist(groupings),,drop=FALSE]
  xPos<-structure(1:nrow(odMat)+cumsum(newGroup)*groupSpace,.Names=rownames(odMat))
  yPos<-structure(ncol(odMat):1,.Names=colnames(odMat))
  odMat[odMat<filterLess]<-filterLess
  if(!is.null(filterMore)){
    anyGreater<-any(odMat>filterMore&!is.na(odMat))
    odMat[odMat>filterMore]<-filterMore
  }
  breaks<-c(-1,seq(log10(filterLess+.5),log10(max(odMat,na.rm=TRUE)+.5),length.out=201))
  cols<-c('white',rev(grDevices::heat.colors(230)[-201:-230]))
  #image(1:nrow(odMat),1:ncol(odMat),log10(odMat),col=cols,breaks=breaks,xaxt='n',yaxt='n',xlab='',ylab='',...)
  graphics::plot(1,1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',xlim=c(.5,max(xPos)+.5),ylim=c(.5,max(yPos)+.5),xaxs='i',yaxs='i',bty='n',...)
  indivCols<-ifelse(is.na(unlist(odMat)),'lightgrey',cols[as.numeric(cut(log10(unlist(odMat)),breaks=breaks))])
  graphics::rect(rep(xPos,ncol(odMat))-.5,rep(yPos,each=nrow(odMat))-.5,rep(xPos,ncol(odMat))+.5,rep(yPos,each=nrow(odMat))+.5,border=NA,xpd=NA,col=indivCols)
  graphics::rect(tapply(xPos,group,min)-.5,graphics::par('usr')[3],tapply(xPos,group,max)+.5,graphics::par('usr')[4],xpd=NA)
  graphics::segments(xPos[-1]-.5,graphics::par('usr')[3],xPos[-1]-.5,graphics::par('usr')[4],col='#00000033')
  graphics::segments(rep(tapply(xPos,group,min)-.5,each=length(yPos)),rep(yPos,3)-.5,rep(tapply(xPos,group,max)+.5,each=length(yPos)),rep(yPos,3)-.5,col='#00000033')
  dnar::multiYAxis(yPos,yLabs,space=.7,axisArgs=list(tcl=-.2),start=.4)
  #axis(1,1:nrow(odMat),rownames(odMat),las=1)
  dnar::slantAxis(1,xPos,rownames(odMat),axisArgs=list(lwd=NA,lwd.ticks=1,tcl=-.3),location=.6)
  #box()
  labPos<-unique(c(200,1000,10000))
  labs<-labPos
  if(!is.null(filterMore)){
    labPos<-unique(c(labPos,filterMore))
    labs<-labPos
    labPos<-labPos[labPos<=filterMore]
    if(anyGreater)labs[labs==filterMore]<-sprintf('>%d',filterMore)
  }
  dnar::insetScale(breaks[-1],cols[-1],c(0.025, 0.035, 0.04, 0.33),at=log10(labPos),labels=labs,main=scaleMain,)
  return(list('x'=xPos,'y'=yPos))
  #graphics::abline(h=2:ncol(odMat)-.5,v=2:nrow(odMat)-.5,col='#00000033')
}


#' Example antibody data
#'
#' A dataset containing raw data read in from the .csv results from an example antibody experiment
#'
#' @format A data frame with 256 rows and 14 cols where the first column gives the patient ID (blanks are filled down), the second column gives the antigen and columns 3-14 give OD values from 96 well plates.
#' @source system.file("data-raw", "makeExample.R", package = "antir")
"antibodyData"
