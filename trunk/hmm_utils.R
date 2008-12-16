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

print_join<-function(nclusts, dir, abnorm=0.9, norm=0.75){
  n<-nclusts%/%2;
  z<-vector(length=nclusts);
  normX<-(1-abnorm)*norm;
#  z[1:n]<-(1-(abnorm+normX))/n;

  z[1:n]<-(1-(abnorm+normX))/(n);
  z[(n+1)]<-normX;
  
  #z[(n+2):nclusts]<-abnorm/n;

  #
  z[(n+2):(nclusts-1)]<-abnorm/(n-1);
  z[nclusts]<-0;
  #
  
  if(dir<0){
    z<-rev(z);
  }  
  cat("[");
  for(i in seq(nclusts-1)){
     cat(z[i],",",sep='');
   }
  cat(z[i+1],sep='');
  cat("],\n");
}

print_abnorm<-function(nclusts, dir, abnorm=0.97){
  n<-nclusts%/%2;
  if(n==1){
    abnorm<-1;
  }
  z<-vector(length=nclusts);
  z[1:(n+1)]<-0;
  z[(n+2):(nclusts-1)]<-getEmisAb(n-1, abnorm);
  z[nclusts]<-abnorm;
  if(dir<0){
    z<-rev(z);
  }  
  cat("[");
  for(i in seq(nclusts-1)){
     cat(z[i],",",sep='');
   }
  cat(z[i+1],sep='');
  cat("],\n");
}

getEmisAb<-function(n, abnorm){
  z<-vector(length=n);
  factRemain<-1-abnorm;
  i<-0;
  if(n>1){
    for(i in seq(n-1)){
      #cat(paste("FactR: ",factRemain,"\n",sep=''));
      z[i]<-factRemain * abnorm;
      factRemain <- factRemain - z[i];
    }
  }
  z[i+1]<-factRemain;
  rev(z);
}


  
print_norm<-function(nclusts, norm=0.6){
  cat("[");
  norm_i<-nclusts%/%2+1;
  neg<-getEmis(nclusts, norm,-1);
  pos<-getEmis(nclusts,norm,1);

 # j<-1;
  for(i in seq(length(neg))){
    cat(neg[i],",",sep='');
  #  j<-j+1;
  }
  cat(norm,",",sep='');
  i<-0;
  if(length(pos)>1){
    for(i in seq(length(pos)-1)){
      cat(pos[i],",",sep='');
    }
  }
  cat(pos[i+1],sep='');
  cat("],\n");
}

getEmis<-function(nclusts, fact, dir){
  factRemain<-(1-fact)/2;
  nRemain<-nclusts%/%2
  z<-vector(length=nRemain);
  i<-0;
  if(nRemain>1){
    for(i in seq(nRemain-1)){
      #cat(paste("FactR: ",factRemain,"\n",sep=''));
      z[i]<-factRemain * fact;
      factRemain <- factRemain - z[i];
    }
  }
  i<-i+1;
  z[i]<-factRemain;
  if(dir<0){
    z<-rev(z);
  }
  z;
    
}


write_emis<-function(nclusts,nstays_pure,norm, abnorm, joinNorm, joinAbnorm){
  #normal
  cat("B = [", sep='');
  print_norm(nclusts,norm);
  #positive
  for(j in seq(nstays_pure)){
    print_abnorm(nclusts, 1, abnorm);
  }
  #joiner pos
  print_join(nclusts, 1, joinAbnorm, joinNorm);
  
  #negative
  for(j in seq(nstays_pure)){
    print_abnorm(nclusts, -1, abnorm);
  }

  #joiner neg
  print_join(nclusts, -1, joinAbnorm, joinNorm);
  
  cat("];\n");
}

print_norm_t<-function(nstay_pure){
  norm_norm<-.999;
  norm_pos<-.0005;

#  norm_norm<-.99999;
#  norm_pos<-.000005;
  cat("[",norm_norm,",",sep='');
  #pos
  for(cj in seq(nstay_pure)){
    if(cj==1){
      cat(norm_pos,",",sep='');
    }
    else{
      cat(0,",",sep='');
    }
  }

  #joiner
  cat(0,",",sep='');
  
  # neg
  for(cj in seq(nstay_pure)){
    if(cj==1){
      cat(norm_pos,",",sep='');
    }
    else{
      cat(0,",",sep='');
    }
  }
  #joiner
  cat(0,",",sep='');

  cat("],\n");
}






print_lr_t<-function(nstates,self,other,joiner,seed=FALSE, lt=data.frame(ab_self=0.84, ab_norm=0.14, ab_joiner=0.01)){
  ab_other<-1-sum(lt$ab_self,lt$ab_norm, lt$ab_joiner);
#  ab_other<-0.01;
  cat("[");
  for(i in seq(nstates)){
    if(i==1){
      cat(lt$ab_norm,",",sep='');
    }
    else{
      if(i==self){
        cat(lt$ab_self,",",sep='');
      }
      else{
        if(i==other){
          cat(ab_other,",",sep='');
        }
        else{
          if(i==joiner){
            cat(lt$ab_joiner,",",sep='');
          }
          else{
            cat(0,",",sep='');
          }
        }
      }
    }
  }
  cat("],\n",sep='');
}


print_abnorm_t<-function(str,nstates,nstay){
  
  for(i in seq(nstay-1)){
 #   cat(str,"\n");
    cat("[");
    for(j in seq(1:str)){
      cat(0,",",sep='');
    }
    cat(1,",",sep='');
    if(str+2 <= nstates){
      for(j in (str+2):nstates){
        cat(0,",",sep='');
      }
    }
    cat("],\n",sep='');
    str<-str+1;
  }
}


print_joiner_t<-function(nstates, self, same, seed=FALSE, jt=data.frame(j_self=0.4, j_same=0.6)){
  cat("[");
  for(i in seq(nstates)){
    if(i==self){
      if(seed){
        cat(0,",",sep='');
      }
      else{
        cat(jt$j_self,",",sep='');
      }
    }
    else{
      if(i==same){
        if(seed){
          cat(1,",",sep='');
        }
        else{
          cat(jt$j_same,",",sep='');
        }
      }
      else{
        cat(0,",",sep='');
      }
    }
  }
  cat("],\n",sep='');
}

write_trans<-function(nstay_pure,lt,jt,seed=FALSE){
#  nstates<-nstay_pure*2+1;
  nstates<-nstay_pure*2+3;
  coords<-matrix(nrow=2,ncol=2);
  coords[1,]<-c(2,2+nstay_pure-1);
  #coords[2,]<-c(coords[1,2]+1,coords[1,2]+nstay_pure);
  coords[2,]<-c(coords[1,2]+2,coords[1,2]+nstay_pure+1);
  
  cat("A = [", sep='');
  print_norm_t(nstay_pure);
  
  print_abnorm_t(coords[1,1],nstates,nstay_pure);
   print_lr_t(nstates=nstates,self=coords[1,2],other=coords[2,1], joiner=coords[2,1]-1,seed=seed, lt=lt);
  #print_joiner_t(nstates, self=coords[1,2]+1, same=coords[1,1], seed=seed);
   print_joiner_t(nstates, self=coords[1,2]+1, same=coords[1,2], seed=seed, jt=jt);
  
  print_abnorm_t(coords[2,1],nstates,nstay_pure);
  print_lr_t(nstates=nstates,self=coords[2,2],other=coords[1,1], joiner=coords[2,2]+1,seed=seed, lt=lt);
  #print_joiner_t(nstates, self=coords[2,2]+1, same=coords[2,1],seed=seed);
  print_joiner_t(nstates, self=coords[2,2]+1, same=coords[2,2],seed=seed, jt=jt);

  cat("];\n", sep='');
}
  
write_dict<-function(nstay_pure){
  cat("dict = {0:0,\n",sep='');
  x<-1;
  i<-1;
  
  for(j in seq(nstay_pure+1)){
#  for(j in seq(nstay_pure)){
    cat(x,":",i,",\n",sep='');
    x<-x+1;
  }
 # cat(x,":",3,",\n",sep='');
 # x<-x+1;
  
  i<-2;
  for(j in seq(nstay_pure+1)){
  #for(j in seq(nstay_pure)){
    cat(x,":",i,",\n",sep='');
    x<-x+1;
  }

  #cat(x,":",4,",\n",sep='');
  #x<-x+1;
  
  cat("};\n",sep='');
}

write_pi<-function(nstay_pure){
  nstates<-nstay_pure*2+3;
  gain_i<-2;
  loss_i<-gain_i+nstay_pure+1;
  cat("pi=[.99,",sep='');
  for(i in 2:(nstates-1)){
    if(i==gain_i || i == loss_i){
      cat(".005,");
    }
    else{
      cat("0,");
    }
  }
  cat("0");
  cat("];\n",sep='');
}
