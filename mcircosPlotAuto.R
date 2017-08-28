#!/usr/bin/env Rscript

###
# Environment
NM<-Sys.getenv(c("NM"))
SRC<-Sys.getenv(c("SRC"))
LENCIRC<-as.numeric(Sys.getenv(c("LENCIRC")))
GENOMEPATH<-Sys.getenv(c("GENOMEPATH"))



# Additional annotation
chromLen<-LENCIRC
annoOther<-read.table(paste(GENOMEPATH,"/AnnoOther",sep=""),header=F)
oril<-annoOther[which(annoOther$V1=="oril"),2]
orih<-annoOther[which(annoOther$V1=="orih"),2]
outerUnit<-1000



rnm<-read.table(paste(NM,".chrM.toLoad",sep=""))
runMchrMf<-rnm[which(rnm$V1=="chrM" & rnm$V3=="+"),]
runMchrMr<-rnm[which(rnm$V1=="chrM" & rnm$V3=="-"),]
totalSignal<-sum(runMchrMf$V4)+sum(runMchrMr$V4)
runMchrMf$V4=(runMchrMf$V4/totalSignal)*1e6
runMchrMr$V4=(runMchrMr$V4/totalSignal)*1e6
NMname<-NM
# Name of eventual output file
NMsvg<-"mCircos.svg"
source(file=paste(SRC,"/pallet.R",sep=""))
source(file=paste(SRC,"/mcircos.R",sep=""))
source(file=paste(SRC,"/chrMplotterMcircosCaller.R",sep=""))
