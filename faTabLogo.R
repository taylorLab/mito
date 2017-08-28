#!/usr/bin/env Rscript

###
# Environment
NM<-Sys.getenv(c("NM"))
SRC<-Sys.getenv(c("SRC"))
LENCIRC<-Sys.getenv(c("LENCIRC"))

###
# Data loading
args = commandArgs(trailingOnly=TRUE)
fname<-args[1]
ti<-args[2]
bg<-args[3]
multiplier<-args[4]
coordStart<-args[5]
fatab<-read.table(file=fname)
seqLen<-nchar(as.character(fatab$V2[1]))
source("../../scripts/pallet.R")
pallet<-c(lineA,lineC,lineG,lineT)

###
# Background nucleotide frequencies
backgroundCounts<-data.frame(equal=rep(10,4))
row.names(backgroundCounts)<-c("A","C","G","T")
backgroundCounts$mm9ChrMfor<-c(5628,3976,2012,4681)
backgroundCounts$mm9ChrMrev<-c(4681,2012,3976,5628)
backgroundCounts$mm9ChrMboth<-backgroundCounts$mm9ChrMfor+backgroundCounts$mm9ChrMrev
backgroundCounts$mm9Nuclear<-c(1335976010,954477394,954477394,1335976010)
backgroundCounts$hg38ChrMfor<-c(5124,5181,2169,4094)
backgroundCounts$hg38ChrMrev<-c(4094,2169,5181,5124)
backgroundCounts$hg38ChrMboth<-backgroundCounts$hg38ChrMfor+backgroundCounts$hg38ChrMrev
backgroundCounts$hg38Nuclear<-c(1560229916,1069882640,1069882640,1560229916)
backgroundRates<-data.frame(apply(backgroundCounts,2,function(x) x/sum(x)))
if(length(which(names(backgroundRates)==bg))<1){
 bg<-"equal"
}

###

###
# Functions
#source("../../scripts/pwmLogo.R")
source(paste(SRC,"/pwmLogo.R",sep=""))

psrIt<-function(psr,multiplier=0,bg="equal",...){
 mmsd<-dim(psr)
 seqLen<-nchar(as.character(psr$V2[1]))
 aah<-matrix(unlist(strsplit(as.character(psr$V2[1:mmsd[1]]),"")),mmsd[1],seqLen,byrow=T)
 nucCountsMatrix<-matrix(rep(0,seqLen*4),4,seqLen)
 nucs<-c("A","C","G","T")
 if(multiplier==1){
  nucCountsMatrix[1,]=apply(aah,2,function(x) sum(psr[which(factor(x,nucs)=="A"),1]))
  nucCountsMatrix[2,]=apply(aah,2,function(x) sum(psr[which(factor(x,nucs)=="C"),1]))
  nucCountsMatrix[3,]=apply(aah,2,function(x) sum(psr[which(factor(x,nucs)=="G"),1]))
  nucCountsMatrix[4,]=apply(aah,2,function(x) sum(psr[which(factor(x,nucs)=="T"),1]))
 }
 else{
  nucCountsMatrix<-apply(aah,2,function(x) table(factor(x,nucs)))
 }
 numFreqMatrix<-apply(nucCountsMatrix,2,function(x) x/sum(x))
 pwmLogo(numFreqMatrix,axis=F,add=T,background=backgroundRates[[bg]],...)
}
###

####
# Plotting
pdf(file=paste(ti,".pdf",sep=""),width=3.5,height=3)

par(mar=c(4,4,4,1))
par(xpd=T)
plot(c(0,seqLen+1),c(0,2),axes=FALSE,type="n",ylab="Bit score",main=ti,xlab="",ylim=c(0,1))
axis(1,tick=T,at=c(1:seqLen),labels=F,lwd.ticks=0.6,tck=-0.05)
mid<-floor(seqLen/2)+1
###
# x-axis units set here
if(coordStart!=1){
 axis(1,tick=T,at=c(1,floor(seqLen/2)+1,seqLen),lwd.ticks=1,labels=c(1-mid,0,mid-1))
} else{
 axis(1,tick=T,at=c(1,floor(seqLen/2)+1,seqLen),lwd.ticks=1)
}
text(0,-0.1,"5'",col="grey",cex=1)
text(seqLen+1,-0.1,"3'",col="grey",cex=1)
axis(2)
mtext(paste("Background:",bg,sep=" "),side=3,line=0.5,col="grey",cex=0.6)
psrIt(fatab,bg=bg,pallet=pallet,multiplier=multiplier)

dev.off()
