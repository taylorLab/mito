#!/usr/bin/env Rscript

###
# Environment
NM<-Sys.getenv(c("NM"))
SRC<-Sys.getenv(c("SRC"))
LENCIRC<-Sys.getenv(c("LENCIRC"))
GENOSHORT<-Sys.getenv(c("GENOSHORT"))
GENOMEPATH<-Sys.getenv(c("GENOMEPATH"))
###
# Data loading
genomeCount<-read.table(paste(GENOMEPATH,"/",GENOSHORT,"WholeGenomeMap100.compo",sep=""),header=T)
mitochCount<-read.table(paste(GENOMEPATH,"/chrM.compo",sep=""),header=T)
###


###
# Build triplets
runs<-c(NM)
triRateVec<-data.frame(runName=runs)
nt<-c("A","C","G","T")
rc<-c("T","G","C","A")
trip<-rep("NNN",64)
tripR<-rep("NNN",64)
co=0
for(i in c(1:4)){
 for(j in c(1:4)){
  for(k in c(1:4)){
    co=co+1
    trip[co]=paste(nt[j],nt[i],nt[k],sep="")
    tripR[co]=paste(rc[k],rc[i],rc[j],sep="")
    triRateVec[trip[co]]=0
  }
 }
}

###

###
# Prepair data frames
monoNucOutput<-data.frame(name=runs,NrA=0,NrC=0,NrG=0,NrT=0,MrA=0,MrC=0,MrG=0,MrT=0,meanNrA=0,meanNrC=0,meanNrG=0,meanNrT=0,meanMrA=0,meanMrC=0,meanMrG=0,meanMrT=0,fractionChrM=0)
nnTemplate<-data.frame(num=c(1:64),nuc=trip,count=0,signal=0,median=0,mean=0)
nrTemplate<-data.frame(num=c(1:64),nuc=tripR,count=0,signal=0,median=0,mean=0)
monoNTempate<-data.frame(num=c(1:4),nuc=nt,count=0,signal=0,median=0,mean=0)
source(file=paste(SRC,"/pallet.R",sep=""))
triRateVecMf<-triRateVec
triRateVecMr<-triRateVec
###


###
# Process library data
for(j in 1:length(runs)){
 nn<-nnTemplate
 #nmF<-monoNTemplateF
 #nmR<-monoNTemplateR
 NM<-runs[j]
 rnm<-read.table(paste(NM,".toLoad",sep=""))
 runMchrMf<-rnm[which(rnm$V1=="chrM" & rnm$V3=="+"),]
 runMchrMr<-rnm[which(rnm$V1=="chrM" & rnm$V3=="-"),]
 rnmRem<-rnm[which(rnm$V1!="chrM"),]
 #for(i in c(1:dim(nmF)[1])){
  #genomeCountConverter<-((i*16))+c(0:15)
  #sites<-which(rnmRem$V5 %in% nn$nuc[genomeCountConverter])
  #monoNTemplate$count[i]=as.numeric(genomeCount[trip[genomeCountConverter]])+as.numeric(genomeCount[tripR[genomeCountConverter]])

  #monoNTemplate$signal[i]=sum()
  #monoNTemplate$median[i]=median()
 for(i in c(1:dim(nn)[1])){
  sites<-which(rnmRem$V5 %in% nn$nuc[i])
  nn$count[i]=as.numeric(genomeCount[trip[i]])+as.numeric(genomeCount[tripR[i]])
  nn$signal[i]=sum(rnmRem$V4[sites])
  nn$median[i]=median(rnmRem$V4[sites])
  nn$mean[i]=nn$signal[i]/nn$count[i]
  triRateVec[j,trip[i]]=nn$mean[i]
  sites<-which(runMchrMf$V5 %in% nn$nuc[i])
  triRateVecMf[j,trip[i]]=sum(runMchrMf$V4[sites])/as.numeric(mitochCount[trip[i]])
  sites<-which(runMchrMr$V5 %in% nn$nuc[i])
  triRateVecMr[j,trip[i]]=sum(runMchrMr$V4[sites])/as.numeric(mitochCount[tripR[i]])
 }
}

###
# Write output files for subsequent collation and comparison
write.table(triRateVec,file="triRateVecGenome.Rtab",quote=F,sep="\t")
write.table(triRateVecMf,file="triRateVecChrMf.Rtab",quote=F,sep="\t")
write.table(triRateVecMr,file="triRateVecChrMr.Rtab",quote=F,sep="\t")
###

# triRateVec<-read.table(file="triRateVecGenome.Rtab",header=T)
# triRateVecMf<-read.table(file="triRateVecChrMf.Rtab",header=T)
# triRateVecMr<-read.table(file="triRateVecChrMr.Rtab",header=T)

###
# Plot the normalised rate vectors
tRVGp<-triRateVec
tRVMfp<-triRateVecMf
tRVMrp<-triRateVecMr
for(i in c(1:dim(tRVGp)[1])){
 tRVGp[i,c(2:65)]=(triRateVec[i,c(2:65)]/sum(triRateVec[i,c(2:65)]))*100
 tRVMfp[i,c(2:65)]=(triRateVecMf[i,c(2:65)]/sum(triRateVecMf[i,c(2:65)]))*100
 tRVMrp[i,c(2:65)]=(triRateVecMr[i,c(2:65)]/sum(triRateVecMr[i,c(2:65)]))*100
}
for(i in c(1:dim(tRVGp)[1])){
 NMl<-as.character(tRVMfp[i,1])
 #pdf(file=paste(NMl,"_trinuc.pdf",sep=""),width=9,height=3)
 pdf(file="trinuc.pdf",width=9,height=3)
 layout(matrix(c(1,2,3),1,3))
 par(mar=c(3,3,2,2))
 yBounds<-c(min(tRVMfp[i,c(1:64)+1],tRVMrp[i,c(1:64)+1],tRVGp[i,c(1:64)+1]),max(tRVMfp[i,c(1:64)+1],tRVMrp[i,c(1:64)+1],tRVGp[i,c(1:64)+1]))
 plot(c(1:64),tRVMfp[i,c(1:64)+1],ylab="Mean normalised signal",xlab="",axes=F,type="n",main=paste(NM,"ChrM",sep=" "),ylim=yBounds)
 shape<-2
 size<-0.5
 points(c(1:16),tRVMfp[i,c(1:16)+1],col=lineA,cex=size,pch=shape)
 points(c(17:32),tRVMfp[i,c(17:32)+1],col=lineC,cex=size,pch=shape)
 points(c(33:48),tRVMfp[i,c(33:48)+1],col=lineG,cex=size,pch=shape)
 points(c(49:64),tRVMfp[i,c(49:64)+1],col=lineT,cex=size,pch=shape)
 shape<-6
 points(c(1:16),tRVMrp[i,c(1:16)+1],col=lineA,cex=size,pch=shape)
 points(c(17:32),tRVMrp[i,c(17:32)+1],col=lineC,cex=size,pch=shape)
 points(c(33:48),tRVMrp[i,c(33:48)+1],col=lineG,cex=size,pch=shape)
 points(c(49:64),tRVMrp[i,c(49:64)+1],col=lineT,cex=size,pch=shape)
 axis(2)
 axis(1,at=c(0.5,16.5,32.5,48.5,64.5),labels=c("","","","",""))
 axis(1,at=c(8,24,40,56),labels=c("A","C","G","T"),tick=F)
 plot(c(1:64),tRVGp[i,c(1:64)+1],ylab="Mean normalised signal",xlab="",axes=F,type="n",main=paste(NM,"Nuclear",sep=" "),ylim=yBounds)
 shape<-21
 points(c(1:16),tRVGp[i,c(1:16)+1],col=lineA,cex=size,pch=shape)
 points(c(17:32),tRVGp[i,c(17:32)+1],col=lineC,cex=size,pch=shape)
 points(c(33:48),tRVGp[i,c(33:48)+1],col=lineG,cex=size,pch=shape)
 points(c(49:64),tRVGp[i,c(49:64)+1],col=lineT,cex=size,pch=shape)
 axis(2)
 axis(1,at=c(0.5,16.5,32.5,48.5,64.5),labels=c("","","","",""))
 axis(1,at=c(8,24,40,56),labels=c("A","C","G","T"),tick=F)
 par(mar=c(4,4,1,1))
 plot(as.numeric(tRVGp[i,c(1:64)+1]),as.numeric(tRVMfp[i,c(1:64)+1]),type="n",axes=F,xlim=yBounds,ylim=yBounds,xlab="Nuclear (normalised signal)",ylab="Mitochondrial (normalised signal)")
 points(tRVGp[i,c(1:16)+1],(tRVMfp[i,c(1:16)+1]+tRVMrp[i,c(1:16)+1])/2,col=lineA,cex=size,pch=shape)
 points(tRVGp[i,c(17:32)+1],(tRVMfp[i,c(17:32)+1]+tRVMrp[i,c(17:32)+1])/2,col=lineC,cex=size,pch=shape)
 points(tRVGp[i,c(33:48)+1],(tRVMfp[i,c(33:48)+1]+tRVMrp[i,c(33:48)+1])/2,col=lineG,cex=size,pch=shape)
 points(tRVGp[i,c(49:64)+1],(tRVMfp[i,c(49:64)+1]+tRVMrp[i,c(49:64)+1])/2,col=lineT,cex=size,pch=shape)
 corval<-cor.test(as.numeric(tRVGp[i,c(1:64)+1]),as.numeric((tRVMfp[i,c(1:64)+1]+tRVMrp[i,c(1:64)+1])/2))
 txpos<-yBounds[1]+((yBounds[2]-yBounds[1])*.6)
 typos<-yBounds[1]+((yBounds[2]-yBounds[1])*.2)
 text(txpos,typos,paste("cor=",round(corval$estimate,3),sep=""))
 axis(1)
 axis(2)
 dev.off()
}

###
