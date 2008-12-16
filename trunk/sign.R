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

perm_pos<-function(ng){
  rows<-sample(seq(nrow(ng)), nrow(ng));
  ng2<-ng;
  ng2$position<-ng[rows,]$position;
  ng2;
}

getpvals<-function(cscores, rscores){
  pvals<-vector(length=length(cscores));
  quants<-quantile(rscores, seq(0.01,1,by=0.01));
  for(i in seq(length(cscores))){
    # default is one
    pvals[i]<-1;
    for(j in seq(100,1)){
      if(cscores[i]>quants[[j]]){
        cat("j: ",j,"\n");
        cat("quants: ",quants[[j]],"\n");
        pvals[i]<-(1-j/100);
        break;
      }
    }
  }
  pvals;
}



#median
scores3<-function(cnvs, scale,ng, fact=2, flip=TRUE){
  cstr<-ng[1,]$position;
  cstp<-ng[nrow(ng),]$position;
  ans<-matrix(nrow=nrow(cnvs), ncol=3);
  scs<-vector(length=nrow(cnvs));

  for(i in seq(nrow(cnvs))){
    marg<-(cnvs[i,]$stp-cnvs[i,]$str)*scale;
    str<-cnvs[i,]$str - marg;
    stp<-cnvs[i,]$stp + marg;
    med<-median(ng[ng$position>=cnvs[i,]$str & ng$position<=cnvs[i,]$stp,]$nimblegen_signal);
    
    if(str<cstr){
      str<-cstr;
    }
    if(stp>cstp){
      stp<-cstp;
    }
    x1<-ng[ng$position>=str & ng$position<cnvs[i,]$str,];
    x2<-ng[ng$position<=stp & ng$position>cnvs[i,]$stp,];
    x<-rbind(x1,x2);
       
    vals2<-x[sign(x$nimblegen_signal)==sign(cnvs[i,]$mean_sig),];
    if(nrow(vals2)<10){
      ans[i,]<-score3(cnvs[i,], med=med, vals=x, fact=fact, flip=FALSE);
    }
    else{
      ans[i,]<-score3(cnvs[i,], med=med, vals=x, fact=fact, flip=flip);
    }
  }
  ans;
}

score3<-function(cnv,med, vals, fact=4, flip=TRUE){
  if(nrow(vals)<2){
    ssdd<-0;
    x<-med;
  }
  else{
    x<-med;
    if(flip){
        vals2<-vals[sign(vals$nimblegen_signal)==sign(cnv$mean_sig),];
        ssdd<-sd(vals2$nimblegen_signal,na.rm=TRUE);
    }
    else{
      ssdd<-sd(vals$nimblegen_signal,na.rm=TRUE);
    }
  }
  sco<-abs(x)*log(cnv$num_probes)-(fact*ssdd);
  c(sco, ssdd, med);
}
