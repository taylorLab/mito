##
# needs..
# NM (run name)
# chromLen
# rnm (main toLoad table sucked in)

svg(file=NMsvg,width=12,height=12)
plot(c(1,outerUnit),c(1,outerUnit),type="n",axes=F,xlab="",ylab="",main="")
xc<-outerUnit/2
yc<-outerUnit/2
yUB<-max(runMchrMr$V4,runMchrMf$V4)
yLB <-0
innerBase<-290
outerBase<-320
yAxisUnitSpan<-180
numGuideLines<-4
ySpan<-abs(yUB-yLB)
yScaleFactor<-ySpan/yAxisUnitSpan
hLineUnits<-floor((ySpan/numGuideLines)/10)*10
scaledLineUnits<-hLineUnits/yScaleFactor
for(o in c(1:4)){
  # Guide concentric rings
  mst.arc(xc=xc,yc=yc,r=innerBase-(o*scaledLineUnits),from=1,to=chromLen,genLength=chromLen,lwd=1,lty=3,col="lightblue")
  mst.arc(xc=xc,yc=yc,r=outerBase+(o*scaledLineUnits),from=1,to=chromLen,genLength=chromLen,lwd=1,lty=3,col="lightblue")
  tx<-mPoint(xc=xc,yc=yc,r=innerBase-(o*scaledLineUnits),point=chromLen*.875,genLength=chromLen)
  points(tx$x[1],tx$y[1],pch=20,cex=3,lwd=2,col="white")
  text(tx$x[1],tx$y[1],hLineUnits*o,col="lightblue",cex=0.6)
  tx<-mPoint(xc=xc,yc=yc,r=outerBase+(o*scaledLineUnits),point=chromLen*.875,genLength=chromLen)
  points(tx$x[1],tx$y[1],pch=20,cex=3,lwd=2,col="white")
  text(tx$x[1],tx$y[1],hLineUnits*o,col="lightblue",cex=0.6)
}
distxyc<-outerBase+((numGuideLines+0.5)*scaledLineUnits)
# Guide crosshairs
points(c(xc,xc),c(yc+distxyc,yc-distxyc),type="l",col="lightblue",lty=3,lwd=1)
points(c(xc+distxyc,xc-distxyc),c(yc,yc),type="l",col="lightblue",lty=3,lwd=1)


# A G C T <- plot order
triRateVec<-data.frame(runName=c(NM))
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

As<-trip[1:16]
Cs<-trip[17:32]
Gs<-trip[33:48]
Ts<-trip[49:64]
AsIndex<-which(runMchrMr$V5 %in% As)
CsIndex<-which(runMchrMr$V5 %in% Cs)
GsIndex<-which(runMchrMr$V5 %in% Gs)
TsIndex<-which(runMchrMr$V5 %in% Ts)
rLB<-innerBase
rUB<-innerBase-yAxisUnitSpan
mitoCircos(xVals=runMchrMr$V2[AsIndex],yVals=runMchrMr$V4[AsIndex],rLowerBound=rLB,rUpperBound=rUB,genLength=chromLen,yUpperBound=yUB,col=lineAa)
mitoCircos(xVals=runMchrMr$V2[GsIndex],yVals=runMchrMr$V4[GsIndex],rLowerBound=rLB,rUpperBound=rUB,genLength=chromLen,yUpperBound=yUB,col=lineGa)
mitoCircos(xVals=runMchrMr$V2[CsIndex],yVals=runMchrMr$V4[CsIndex],rLowerBound=rLB,rUpperBound=rUB,genLength=chromLen,yUpperBound=yUB,col=lineCa)
mitoCircos(xVals=runMchrMr$V2[TsIndex],yVals=runMchrMr$V4[TsIndex],rLowerBound=rLB,rUpperBound=rUB,genLength=chromLen,yUpperBound=yUB,col=lineTa)

AsIndex<-which(runMchrMf$V5 %in% As)
CsIndex<-which(runMchrMf$V5 %in% Cs)
GsIndex<-which(runMchrMf$V5 %in% Gs)
TsIndex<-which(runMchrMf$V5 %in% Ts)
rLB<-outerBase
rUB<-outerBase+yAxisUnitSpan
mitoCircos(xVals=runMchrMf$V2[AsIndex],yVals=runMchrMf$V4[AsIndex],rLowerBound=rLB,rUpperBound=rUB,genLength=chromLen,yUpperBound=yUB,col=lineAa)
mitoCircos(xVals=runMchrMf$V2[GsIndex],yVals=runMchrMf$V4[GsIndex],rLowerBound=rLB,rUpperBound=rUB,genLength=chromLen,yUpperBound=yUB,col=lineGa)
mitoCircos(xVals=runMchrMf$V2[CsIndex],yVals=runMchrMf$V4[CsIndex],rLowerBound=rLB,rUpperBound=rUB,genLength=chromLen,yUpperBound=yUB,col=lineCa)
mitoCircos(xVals=runMchrMf$V2[TsIndex],yVals=runMchrMf$V4[TsIndex],rLowerBound=rLB,rUpperBound=rUB,genLength=chromLen,yUpperBound=yUB,col=lineTa)


##
# Chromosome annotation
mst.arcPoly(xc=xc,yc=yc,r=310,from=1,to=chromLen,length=chromLen,width=4,col="lightgrey")
mst.arcPoly(xc=xc,yc=yc,r=300,from=1,to=chromLen,length=chromLen,width=4,col="lightgrey")
#chrAnnoCds<-read.table(file="NC_005089.cds",header=F)
#chrAnnoTrn<-read.table(file="NC_005089.transcript",header=F)
chrAnnoCds<-read.table(file=paste(GENOMEPATH,"/AnnoCds",sep=""),header=F)
chrAnnoTrn<-read.table(file=paste(GENOMEPATH,"/AnnoTrn",sep=""),header=F)


ppos<-((chrAnnoCds$V3-chrAnnoCds$V2)/2)+chrAnnoCds$V2
tx<-mPoint(xc=xc,yc=yc,r=470,point=ppos,genLength=chromLen)
for(g in c(1:dim(chrAnnoCds)[1])){
  caRadius<-310
  if(chrAnnoCds$V5[g] == "-"){
      caRadius=300
  }
  mst.arcPoly(xc=xc,yc=yc,r=caRadius,from=chrAnnoCds$V2[g],to=chrAnnoCds$V3[g],width=10,length=chromLen,col="darkgrey")
  text(tx$x[g],tx$y[g],chrAnnoCds$V4[g],col="darkgrey",cex=0.8)
}
for(g in c(1:dim(chrAnnoTrn)[1])){
  caRadius<-310
  if(chrAnnoTrn$V5[g] == "-"){
    caRadius=300
  }
  mst.arcPoly(xc=xc,yc=yc,r=caRadius,from=chrAnnoTrn$V2[g],to=chrAnnoTrn$V3[g],width=10,length=chromLen,col="gold")
}
# D-loop
#mst.arcPoly(xc=xc,yc=yc,r=300,from=576,to=16024,width=10,length=chromLen,col="red")
#CSB1
#mst.arcPoly(xc=xc,yc=yc,r=310,from=219,to=233,width=10,length=chromLen,col="blue")
#CSB1
#mst.arcPoly(xc=xc,yc=yc,r=310,from=16087,to=16101,width=10,length=chromLen,col="blue")

#OriH
tx<-mPoint(xc=xc,yc=yc,r=292,point=orih,genLength=chromLen)
ty<-mPoint(xc=xc,yc=yc,r=308,point=orih,genLength=chromLen)
tz<-mPoint(xc=xc,yc=yc,r=300,point=orih-100,genLength=chromLen)
tyt<-mPoint(xc=xc,yc=yc,r=470,point=191,genLength=chromLen)
polygon(c(tx$x[1],ty$x[1],tz$x[1]),c(tx$y[1],ty$y[1],tz$y[1]),col="red",border=F)
text(tyt$x[1],tyt$y[1],expression('Ori'[H]),col="red",cex=0.8)

#OriL
tx<-mPoint(xc=xc,yc=yc,r=302,point=oril,genLength=chromLen)
ty<-mPoint(xc=xc,yc=yc,r=318,point=oril,genLength=chromLen)
tz<-mPoint(xc=xc,yc=yc,r=310,point=oril+100,genLength=chromLen)
tyt<-mPoint(xc=xc,yc=yc,r=470,point=oril,genLength=chromLen)
polygon(c(tx$x[1],ty$x[1],tz$x[1]),c(tx$y[1],ty$y[1],tz$y[1]),col="purple",border=F)
text(tyt$x[1],tyt$y[1],expression('Ori'[L]),col="purple",cex=0.8)


text(10,980,paste(NM,": ",NMname,sep=""),pos=4,col="darkgrey")
legend(50,100,c("rA","rG","rC","rT/U"),col=c(lineAa,lineGa,lineCa,lineTa),pch=20,bty="n",text.col="white",pt.cex=3)
legend(50,100,c("rA","rG","rC","rT/U"),col=c(lineA,lineG,lineC,lineT),pch=20,bty="n",text.col="darkgrey",pt.cex=1.5)

dev.off()
