
################################################################################
#
#       Copyright (C) Patrick Cahan 2007-2009
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

w_plotCNVR<-function(cnvrGenos, cnvrid, marg=5){
  cnvr<-cnvrGenos[['cnvr']] [ cnvrGenos[['cnvr']]$cnvrid==cnvrid,];
  coding<-cnvrGenos[['coding']] [ cnvrGenos[['coding']]$cnvrid==cnvrid,];

  snames<-colnames(coding)[2:ncol(coding)];
  ng<-multiSample(snames,chr=cnvr$chr);

  str<-cnvr$str-(cnvr$length*marg);
  stp<-cnvr$stp+(cnvr$length*marg);

  ng<-sliceNG(ng, data.frame(str=str,stp=stp));

  plotCNVR(ng, cnvr, coding)
}

plotCNVR<-function(ng,cnvr,coding,smooth=FALSE){
  positions<-as.numeric(colnames(ng));
  snames<-colnames(coding)[2:ncol(coding)];
  ylim<-c(min(ng)-.2, max(ng)+.2);
  
  colls<-colors();
  greens<-colls[grep("green", colls)];
  reds<-colls[grep("red", colls)];

  ggg<-1;
  rrr<-1;
  bbb<-1;

  for(i in seq(nrow(ng))){
    sname<-snames[i];
    col1<-"gray";     
    bbb<-bbb+1;
    if(coding[sname]<0){
      col2<-reds[rrr];
      rrr<-rrr+1;
    }
    else{
      if(coding[sname]>0){
         col2<-greens[ggg];
         ggg<-ggg+1;
       }
      else{
        col2<-col1;
      }
    }

    colLeft<-rep(col1, length(positions[positions<cnvr$str]));
    colIn<-rep(col2, length(positions[ (positions>= cnvr$str) &  (positions<=cnvr$stp)]));    
    colRight<-rep(col1, length(positions[positions>cnvr$stp]));
    col2<-c(colLeft, colIn, colRight);
    plot(as.numeric(colnames(ng)),ng[i,],col=col2, pch=16, cex=.5, ylim=ylim,xlab='', ylab='');
    par(new=TRUE);
  }  
  plot_vertical_lines(cnvr, ylim[1], ylim[2]);
  mtext(paste("Chr",cnvr$chr,sep=''), 1, 2.5);
  mtext("Log2(ratio)",2,2.5);
  mtext(paste("CNVR ID:",cnvr$cnvrid), 3,.5)
}
