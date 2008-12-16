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

# return:
#        (1) a table of cnvrs
#        (2) a list of tables, each containing the correpsonding cnvs
#

w_makeCNVRs<-function(calls,strID=1){
  chrs<-unique(as.vector(calls$chr));
  cnvrs<-data.frame();
  cnvBag<-list();
  
  for(chr in chrs){
    callsX<-calls[calls$chr==chr,];
    x<-makeCNVRs(callsX, strID);
    for(cnvrid in x[[1]]$cnvrid){
      cnvBag[[cnvrid]]<-x[[2]][[cnvrid]];
    }
    cnvrs<-rbind(cnvrs,x[[1]]);
    strID<-(max(cnvrs$cnvrid)+1);
  }
  list(cnvrs=cnvrs, cnvrBag=cnvBag);
}

# assume already broken down by chr
makeCNVRs<-function(calls, strID=1){

  cnvBag<-list();
  cur_id<-strID;
  n_added<-0;
  calls<-calls[order(calls$str, calls$stp),];

  if(nrow(calls)>1){
    i<-1;
    while(TRUE){
      # last row
      if(i == nrow(calls)){
        if(n_added==0){
          cnvBag[[cur_id]]<-calls[i,];
        }        
        break;
      }
      else{
        if(ol(calls[i,], calls[(i+1),])){
          callX<-calls[(i+1),];          
          # update cnvrBag
          if(n_added==0){
            cnvBag[[cur_id]]<-data.frame(calls[i:(i+1),]);
            n_added<-2;
          }
          else{
            cnvBag[[cur_id]]<-rbind(cnvBag[[cur_id]], callX);
            n_added<-n_added+1;
          }
          
          # if next to last row
          if(i == (nrow(calls)-1)){
            calls<-calls[1:i,];
          }
          else{
            calls<-calls[c(1:i,(i+2):nrow(calls)),];
          }
          calls[i,]<-mergeCalls(calls[i,], callX);
        }
        # Doesn't overlap next call
        else{
              
          if(n_added==0){
            cnvBag[[cur_id]]<-data.frame(calls[i,]);
          }
          n_added<-0;
          cur_id<-cur_id+1;
          i<-i+1;
        }
      }
    }
  }

  else{
    cnvBag[[cur_id]]<-calls[1,];
  }
  
  cnvrids<-strID:cur_id;
  calls<-cbind(calls, cnvrid=cnvrids);
  ave_score<-vector();
  lens<-vector();
  for(i in seq(nrow(calls))){
    ave_score<-append(ave_score, mean(cnvBag[[calls[i,]$cnvrid]]$score));
    lens<-append(lens, calls[i,]$stp-calls[i,]$str+1);
  }

  calls<-cbind(calls, ave_score=round(ave_score,3));
  calls<-cbind(calls, length=lens);
  calls<-calls[,c('cnvrid','chr','str','stp','length','ave_score')];
  rownames(calls)<-calls$cnvrid;
  list(cnvrs=calls,cnvBag=cnvBag);
  
}


filterCalls<-function(calls,lossScore=1,gainScore=1){
  callsG<-calls[calls$mean_sig>0 & calls$score>gainScore,];
  callsL<-calls[calls$mean_sig<0 & calls$score>lossScore,];
  rbind(callsG, callsL);
}
  

mergeCalls<-function(c1, c2){
  ans<-c1;
  ans$str<-min(c1$str, c2$str);
  ans$stp<-max(c1$stp, c2$stp);
  ans;
}

ol<-function(c1, c2){
  if( (c1$str <= c2$stp) & (c1$stp >= c2$str)){
    ans<-TRUE;
  }
  else{
    ans<-FALSE;
  }
  ans;
}

