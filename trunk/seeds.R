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


#figure out how many probes in each seed
countSeedProbes<-function(ng, seeds){
  nprobes<-vector(length=nrow(seeds));
  for(i in seq(nrow(seeds))){
    seed<-seeds[i,];
    nprobes[i]<-nrow(ng[ng$position>=seed$str & ng$position<=seed$stp,]);
  }
  nprobes;
}

# margin is a factor of the number of probes in the seed
expandSeeds<-function(ng,seeds,marg){
  if(nrow(seeds)>0){
    nprobes<-nrow(ng);
    for(i in seq(nrow(seeds))){
      marg2<-floor(seeds[i,]$num_probes*marg);
      ngStrI<-which(ng$position==seeds[i,]$str);
      if(ngStrI<=marg2){
        seeds[i,]$str<-ng[1,]$position;
      }
      else{
        seeds[i,]$str<-ng[ngStrI-marg2,]$position;
      }

      ngStpI<-which(ng$position==seeds[i,]$stp);
      if( ngStpI > nprobes-marg2 ){
        seeds[i,]$stp<-ng[nprobes,]$position;
      }
      else{
        seeds[i,]$stp<-ng[(ngStpI+marg2),]$position;
      }
    }    
  }  
  seeds;
}


mergeSeeds<-function(seeds){
  i<-1;
  while(i<nrow(seeds)){
    if(i+1>nrow(seeds)){
      break;
    }
    else{
      if(overlap(seeds[i,], seeds[i+1,])){
        seeds<-combSeeds(seeds, i);
        next;
      }
      else{
        i<-i+1;
      }
    }
  }
  seeds;
}

overlap<-function(c1,c2){
  if((c1$str<=c2$stp & c1$stp>=c2$str) ||
     (c1$stp>=c2$str & c1$stp<=c2$stp) ){
    TRUE;
  }
  else{
    FALSE;
  }
}

combSeeds<-function(seeds, i){
  str<-seeds[i,]$str;
  stp<-seeds[i+1,]$stp;
  if(i+1==nrow(seeds)){
    out<-seeds[1:i,];
  }
  else{
    out<-seeds[c(1:i,(i+2):nrow(seeds)),];
  }
  out[i,]$str<-str;
  out[i,]$stp<-stp;
  out;
}

