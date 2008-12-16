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


# compute CNVr genotypes for supplied samples

# Determines the optimal clustering of CNV-region as determined by
# the silhouette and agreement between clustering and segment calls
# returns list of:
#    (1) cnvr: cnvr coordinates, score, etc
#    (2) genotype: clustering assigned labels
#    (3) coding: i.e. -1, 0, 1 for losses, normal, gains, respectively
#    (4) mean_sig: mean log2(ratio) of each sample in cnvr
# Each element can be accessed by usual R naming.So, x[['cnvr']], will return a table of CNVRs.

w_genotype<-function(snames,cnvBag,cnvrs,nclusts=2:7){

  cnvrs<-cnvrs[order(cnvrs$chr, cnvrs$str),];

  cnvrsOut<-data.frame();
  genotypes<-data.frame();
  coding<-data.frame();
  mean_sigs<-data.frame();

  chrs<-unique(as.vector(cnvrs$chr));
  for(chr in chrs){
    x<-genotypeChr(snames,cnvBag,cnvrs,chr,nclusts);
    cnvrsOut<-rbind(cnvrsOut, x$cnvrs);
    genotypes<-rbind(genotypes, x$genotypes);
    coding<-rbind(coding, x$coding);
    mean_sigs<-rbind(mean_sigs, x$mean_sigs);
  }
  
  list(cnvrs=cnvrsOut, genotypes=genotypes, coding=coding, mean_sigs=mean_sigs);
}



#### HELPER fucntions ##### 


# in: ngData, cnvs, nclusts
genotypeCNVR<-function(ngData,cnvr,cnvs,nclusts=2:7){
  array_ids<-as.vector(rownames(ngData));

  if(max(nclusts)> (length(array_ids)-1)){
    nclusts<-2:(length(array_ids)-1);
  }
  
  aveWidths<-vector(length=length(nclusts));
  goodnesses<-vector(length=length(nclusts));
  genos<-list();
  meanSigs<-list();
  i<-1;
  for(nclust in nclusts){
    x<-pam(ngData, nclust);
    labels<-data.frame(array_id=array_ids, label=x$clustering);
    aveWidths[i]<-x$silinfo$avg.width;
    meanSigs[[i]]<-getMeanSigs(ngData, labels);
    refGrp<-detRefGroup(labels$label,meanSigs[[i]]);
    genos[[i]]<-makeCNVrGeno(labels,cnvs,refGrp);
    goodnesses[i]<-agreement(cnvs,genos[[i]],refGrp);
    i<-i+1;
  }

  # select the best fit
  best<-which.max(goodnesses*aveWidths);
  mean_sigs<-meanSigs[[best]];
  genotypes<-getGenotypes(mean_sigs, genos[[best]]$label, array_ids);
  list(ave_width=aveWidths[[best]],
       agreement=goodnesses[[best]],
       nclusts=nclusts[best],
       genotypes=genotypes);
}

# compute CNVr genotypes for supplied samples and chr
genotypeChr<-function(snames,cnvBag,cnvrs,chr,nclusts=2:7){
  cnvrs<-cnvrs[cnvrs$chr==chr,];
  genotypes<-data.frame();
  coding<-data.frame();
  mean_sigs<-data.frame();
  cnvrsOut<-data.frame();
  
  ng<-multiSample(snames,chr);

  for( i in seq(nrow(cnvrs))){
    cnvr<-cnvrs[i,];
    x<-sliceNG(ng,cnvr);
    cnvs<-cnvBag[[cnvr$cnvrid]];
    genos<-genotypeCNVR(x,cnvr,cnvs,nclusts);
    ave_width<-genos$ave_width;
    cnvr<-cbind(cnvr, ave_width=ave_width);
    cnvrsOut<-rbind(cnvrsOut,cnvr);
    genotypes<-rbind(genotypes, data.frame(t(genos$genotypes$genos)));
    coding<-rbind(coding, data.frame(t(genos$genotypes$coding)));
    mean_sigs<-rbind(mean_sigs,data.frame(t(genos$genotypes$meanSigs)));
  }

  genotypes<-cbind(data.frame(cnvrid=cnvrs$cnvrid), genotypes);
  colnames(genotypes)<-c('cnvrid',snames);

  coding<-cbind(data.frame(cnvrid=cnvrs$cnvrid), coding);
  colnames(coding)<-c('cnvrid',snames);

  mean_sigs<-cbind(data.frame(cnvrid=cnvrs$cnvrid), mean_sigs);
  colnames(mean_sigs)<-c('cnvrid',snames);

  list(cnvrs=cnvrsOut, genotypes=genotypes, coding=coding, mean_sigs=mean_sigs);
}

sliceNG<-function(ng, loc){
  thisNG<-ng[ng$position>=loc$str & ng$position<=loc$stp,];
  x<-t(thisNG[,c(2:ncol(thisNG))]);
  colnames(x)<-thisNG$position;
  x;
}

# return the cluster label corresponding to the inds with the log2ratio closest to zero
# assume correct ordering 
detRefGroup<-function(clusterLabels, meanSigs){
  x<-which.min(abs(meanSigs));
  clusterLabels[x];
}

agreement<-function(hcalls, labels, refGrp){
  hcs<-hcalls$name;
  x<-intersect(hcs, labels[labels$label==refGrp,]$array_id);
  1-(length(x)/nrow(hcalls));
}

getMeanSigs<-function(ngData, labels){
  means<-vector(length=nrow(labels));
  for(i in seq(nrow(labels))){
    means[i]<-mean(ngData[i,])
  }
  means;
}
  
# disallows the formation of de novo genotype groups
makeCNVrGeno<-function(labels,hCalls,refGrp){
  grps<-unique(labels$label);
  for(grp in grps){
    if(grp!=refGrp){
      wx<-which(labels$label==grp);
      x<-labels[wx,];
      ol<-length(intersect(x$array_id, hCalls$name));
      if(ol==0){
        labels[wx,]$label<-refGrp;
      }
    }
  }
  labels;
}

# assumes only one 'normal' genotype
getGenotypes<-function(meanSigs, genos, array_ids){
  alleles<-unique(genos);
  df<-data.frame(allele=alleles,
                 gMeans=vector(length=length(alleles)),
                 coding=vector(length=length(alleles)));
  
  for(i in seq(length(alleles))){
    df$gMeans[i]<-mean(meanSigs[which(genos==alleles[i])]);
  }
  # "normal"
  norm<-which.min(abs(df$gMeans));
  df$coding[norm]<-0;
  
  losses<-df[setdiff(which(df$gMeans<0),norm),];
  if(nrow(losses)>0){
    vals<- -1*order(abs(losses$gMeans));
    df[rownames(losses),]$coding<-vals;    
  }

  gains<-df[setdiff(which(df$gMeans>0),norm),];

  if(nrow(gains)>0){
    vals<- order(abs(gains$gMeans));
    df[rownames(gains),]$coding<-vals;    
  }

  ans<-data.frame(array_id=array_ids,
                  meanSigs=meanSigs,
                  genos=genos,
                  coding=vector(length=length(array_ids)));
  for(i in seq(nrow(ans))){
    ans$coding[i]<-df[which(df$allele==genos[i]),]$coding;
  }
  ans;
}
