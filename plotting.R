
################################################################################
#
#       Copyright (C) Patrick Cahan 2007-2008
#
#       Contact: pcahan@wustl.edu
#
#
#       This library is free software; you can redistribute it and/or
#       modify it under the terms of the GNU Library General Public
#       License as published by the Free Software Foundation; either
#       version 2 of the License, or (at your option) any later version.
#
#       This library is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#       Library General Public License for more details.
#
#       You should have received a copy of the GNU Library General Public
#       License along with this library; if not, write to the Free
#       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
#################################################################################

plotReg<-function(sampleName,cnvs=data.frame(),marg=5,ylim=c(-3,3),res=50){

  main<-sampleName;
  chr<-cnvs[1,]$chr
  ng<-readNG(sampleName,chr=chr);

  marg<-(cnvs[1,]$stp-cnvs[1,]$str+1) * marg;

  if(nrow( ng[ng$position >= (cnvs[1,]$str-marg),])>0){
    left<-cnvs[1,]$str-marg;
  }
  else{
    left<-ng[1,]$position;
  }

  if(nrow( ng[ng$position <= (cnvs[1,]$stp+marg),])>0){
    right<-cnvs[1,]$stp+marg;
  }
  else{
    right<-ng[nrow(ng),]$position;
  }

  ng<-ng[ng$position>=left & ng$position<=right,];
    

  # make colors vector
  clrs<-vector();
  cls<-binColors(ng,res=res);
  x<-rgb(cls,cls,1);
  clrs<-append(clrs, x);
  plot(ng$position,ng$nimblegen_signal, main=main,xlab=paste("Chr",chr,sep=''),ylim=ylim,pch='.',cex=2.45,col=clrs, ylab="log2-ratio");
  lines(x=c(0,nrow(ng)), y=c(0,0));
  if(nrow(cnvs)>0){
    plot_cnv_lines(cnvs, cl='red', lty=1);
  }
  
}

binColors<-function(ng, res=100){
  vals<-seq(from=0,to=1.5,length.out=res);
  cls<-seq(from=.5,to=0,length.out=res);
  bins<-vector(length=nrow(ng));
  i<-2;
  while(i < length(cls)){
    bins[which(abs(ng$nimblegen_signal)>=vals[i-1] & abs(ng$nimblegen_signal)<vals[i])]<-cls[i-1];
    i<-i+1;
  }
  bins[which(abs(ng$nimblegen_signal)>=vals[i])]<-cls[i];
  bins;
}

plot_cnv_lines<-function(hmm_cnvs, cl=4, lty=1){
  for(i in seq(nrow(hmm_cnvs))){
    seg<-hmm_cnvs[i,];
    lines(c(seg$str, seg$stp), c(seg$mean_sig, seg$mean_sig), col=cl, lty=lty, lwd=2);
  }
}

plot_vertical_lines<-function(segs, mn, mx){
  for(i in seq(nrow(segs))){
    plot_vertical_line(segs[i,],mn=mn,mx=mx);
  }
}

plot_vertical_line<-function(seg,mn,mx){
  lines(c(seg$str, seg$str), c(mn,mx), col='black', lwd=1, lty=2);
  lines(c(seg$stp, seg$stp), c(mn,mx), col='black', lwd=1,lty=2);
}
