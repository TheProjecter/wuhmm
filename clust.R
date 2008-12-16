
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

# adds a rank (ascending) column to ng based on clustering
addClustRank<-function(ng,num_clusts,diver=FALSE){
  ranks<-vector(length=nrow(ng));
  if(diver){
    similar<-ng[ng$state==1 | ng$state==2,];
    divergent<-ng[ng$state==0,];
    if( (nrow(similar) == nrow(ng)) ||
       (nrow(divergent) == nrow(ng)) ||
       (nrow(similar)<num_clusts) ||
       (nrow(divergent)<num_clusts) ){
      diver<-FALSE;
      ranks<-clustPosNeg(ng,num_clusts);
    }
    else{
      cl_s<-clustPosNeg(similar,num_clusts);
      cl_d<-clustPosNeg(divergent,num_clusts);
      ranks[as.integer(rownames(similar))]<-cl_s;
      ranks[as.integer(rownames(divergent))]<-cl_d;
    }
  }
  else{
    ranks<-clustPosNeg(ng,num_clusts);
  }
  ng<-cbind(ng,ranks);
  #write.table(ng, file=fname, sep="\t", col.names=FALSE, row.names=FALSE,quote=FALSE);
}

# returns ranks of values in tab (ng)
clustPosNeg<-function(tab,nclusts){
  centers<-vector(length=nclusts, mode="integer");
  ranks<-vector(length=nrow(tab),mode="integer");
  
  #divide by mean value
  mn<-mean(tab$nimblegen_signal);
  neg<-tab[tab$nimblegen_signal<mn,];
  indexn<-which(tab$nimblegen_signal<mn);

  pos<-tab[tab$nimblegen_signal>=mn,];
  indexp<-which(tab$nimblegen_signal>=mn);

  #cluster negs
  cln<-clara(neg$nimblegen_signal,nclusts%/%2 + 1);
  
  # set centers
  centers[1:(nclusts%/%2+1)]<-cln$medoids;

  #set cluster
  ranks[indexn]<-translateRank(clus=cln, nclusts=nclusts, dir=-1);
  
  # cluster pos
  clp<-clara(pos$nimblegen_signal,nclusts%/%2 + 1);

  #set cluster
  ranks[indexp]<-translateRank(clus=clp, nclusts=nclusts, dir=1);
  ranks;
}

# translates cluster labels to ranks
translateRank<-function(clus, nclusts, dir){
  center<-(nclusts%/%2)+1;
  ranks<-vector(length=length(clus$clustering), mode="integer");
  # get norm index
  if(dir<0){
    norm_i<-which.max(clus$medoids);
    str<-0;
#    cat(str,"\n");
  }
  else{
    norm_i<-which.min(clus$medoids);
    str<-center;
#    cat(str,"\n");
  }

  ranks[which(clus$clustering==norm_i)]<-center;
  meds<-rem(clus$medoids, clus$medoids[norm_i]);
  medsS<-sort(meds);
  l<-1;
  for(med in clus$medoids){
    if(l!=norm_i){
      rnk<-which(medsS==med) + str;
#      cat(med,"\t");
#      cat(rnk,"\n");
      ranks[which(clus$clustering==l)]<-rnk
    }
    l<-l+1;
  }
  ranks;
}

# removes val from ar
rem<-function(ar, val){
  ans<-vector();
  for(i in ar){
    if(i!=val){
      ans<-append(ans, i);
    }
  }
  ans;
}
