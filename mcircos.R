#
# Martin Taylor Jan 2017
# Simplified version of Circos style plotting

mitoCircos<-function(xVals,yVals,col="#FF000040",rLowerBound,rUpperBound,yLowerBound=0,yUpperBound=NA,genLength=16299){
  rScaledSignalSpan<-rUpperBound-rLowerBound
  if(is.na(yUpperBound)){
    yUpperBound=max(yVals,na.rm=T)
  }
  yBounds<-c(yLowerBound,yUpperBound)
  ySignalSpan<-yBounds[2]-yBounds[1]
  scaledSignal<-(yVals/ySignalSpan)*rScaledSignalSpan
  # The genLength-xVals reverses the order so plots clockwise rather than anticlockwise.
  # *2*pi scales fractions to angles of circle.
  xValMod<-((genLength-xVals)/genLength)*2*pi
  # Moves origin from 3 o'clock to 12 o'clock
  xValMod<-xValMod+(.25*2*pi)
  xAnchorEnd<-xc+cos(xValMod)*rLowerBound
  yAnchorEnd<-yc+sin(xValMod)*rLowerBound
  xRaggedEnd<-xc+cos(xValMod)*(rLowerBound+scaledSignal)
  yRaggedEnd<-yc+sin(xValMod)*(rLowerBound+scaledSignal)
  for(i in c(1:length(scaledSignal))){
    points(c(xAnchorEnd[i],xRaggedEnd[i]),c(yAnchorEnd[i],yRaggedEnd[i]),type="l",col=col,lwd=1,lend=2)
  }
}

mst.arcPoly<-function(xc,yc,r,from,to,length,col="lightblue",width=2,lty=1,lend=1,genLength=16299){
  # xc and yc are the xy for circle midpoint
  whalf<-width/2
  ang.d<-abs(from-to)
  pix.n<-ang.d*5
  if(pix.n<2){
    pix.n<-2
  }
  #ang.seq<-rev(seq(w1,w2,length.out=pix.n))
  ####
  xValMod<-((genLength-c(from,to))/genLength)*2*pi
  # Moves origin from 3 o'clock to 12 o'clock
  xValMod<-xValMod+(.25*2*pi)
  ang.seq<-seq(xValMod[1],xValMod[2],length.out=pix.n)
  ####
  #ang.seq<-ang.seq/360*2*pi
  fan.i.xOuter<-xc+cos(ang.seq)*(r+whalf)
  fan.i.xInner<-xc+cos(ang.seq)*(r-whalf)
  fan.i.yOuter<-yc+sin(ang.seq)*(r+whalf)
  fan.i.yInner<-yc+sin(ang.seq)*(r-whalf)
  #lines(fan.i.xOuter,fan.i.yOuter,col="darkgreen",lwd=lwd,type="l",lend=2,lty=lty)
  #lines(fan.i.xInner,fan.i.yInner,col="green",lwd=lwd,type="l",lend=2,lty=lty)
  polygon(c(fan.i.xOuter,rev(fan.i.xInner)),c(fan.i.yOuter,rev(fan.i.yInner)),col=col,border=NA)
}
#mst.arcPoly(xc=xc,yc=yc,r=500,from=1,to=5000,length=16299,col="#FF000044")


mst.arc<-function(xc,yc,r,from,to,length,col="lightblue",lwd=1,lty=1,lend=1,genLength=16299){
  # xc and yc are the xy for circle midpoint
  ang.d<-abs(from-to)
  pix.n<-ang.d*5
  if(pix.n<2){
    pix.n<-2
  }
  #ang.seq<-rev(seq(w1,w2,length.out=pix.n))
  ####
  xValMod<-((genLength-c(from,to))/genLength)*2*pi
  # Moves origin from 3 o'clock to 12 o'clock
  xValMod<-xValMod+(.25*2*pi)
  ang.seq<-seq(xValMod[1],xValMod[2],length.out=pix.n)
  ####
  #ang.seq<-ang.seq/360*2*pi
  fan.i.x<-xc+cos(ang.seq)*r;
  fan.i.y<-yc+sin(ang.seq)*r
  lines(fan.i.x,fan.i.y,col=col,lwd=lwd,type="l",lend=2,lty=lty)
}
#mst.arc(xc=xc,yc=yc,r=500,from=1,to=5000,length=16299,col="#FF000044")

mPoint<-function(xc,yc,r,point,genLength=16299){
  pointMod<-((genLength-point)/genLength)*2*pi
  pointMod<-pointMod+(.25*2*pi)
  rt<-data.frame(x=rep(0,length(point)),y=rep(0,length(point)))
  rt$x<-xc+cos(pointMod)*r
  rt$y<-yc+sin(pointMod)*r
  return(rt)
}
