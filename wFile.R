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

makeChrs<-function(chrs=c(1:19),x=TRUE, y=TRUE, prepend0=FALSE){
  ans<-vector();

  if(prepend0){
    for(chr in chrs){
      if(chr<10){
        a<-paste('chr0',chr,sep='');
      }
      else{
        a<-paste('chr',chr,sep='');
      }
      ans<-append(ans,a);
    }
  }
  else{
    for(chr in chrs){
      ans<-append(ans,paste('chr',chr,sep=''));
    }
  }
    
  if(x){
    ans<-append(ans, 'chrX');
  }
  if(y){
    ans<-append(ans, 'chrY');
  }  
  ans;
}

# reads a ng data file
# in: sample
# if no chr supplied then get whole genome file
readNG<-function(sname,chr=0, suffix='.wutxt'){

  if(chr==0){
    fname<-paste(sname,"_genome",suffix,sep='');
  }
  else{
    fname<-paste(sname,"_chr",chr,suffix,sep='');
  }
  x<-read.delim(fname,as.is=TRUE);
  colnames(x)<-c('position', 'nimblegen_signal');
  x;
}

colSpec<-function(ncols,chr,position,ng){
  ans<-list();
  for(i in seq(length(ncols))){
    ans[i]<-'NULL';
  }  
  ans[chr]<-'character';
  ans[position]<-'integer';
  ans[ng]<-'numeric';
  ans;
}

stripNGfile<-function(fnameIn,sname,iPos,iSignal,iChr,ncols=14, autos=c(1:19),sex=c(TRUE, TRUE), prepend0=FALSE,suffix='.wutxt'){
  cs<-colSpec(ncols,iChr,iPos, iSignal);
  cat(fnameIn,"\n");
  system.time(xt<-read.table(fnameIn,colClass=cs, skip=1));
  cnames<-data.frame(names=c('position', 'chr', 'ng'), index=c(iPos,iChr, iSignal));
  colnames(xt)<-cnames$names[order(cnames$index)];

  xt$chr<-tolower(xt$chr);
  
  chrs<-makeChrs(chrs=autos,x=sex[1], y=sex[2], prepend0=prepend0);
  for(chr in chrs){
    fnameOut<-paste(sname,"_",chr,suffix,sep='');
    out<-xt[xt$chr==chr,c("position", "ng")];
    out<-out[order(out$position),]
    write.table(out, file=fnameOut,sep="\t",row.names=FALSE, col.names=FALSE, quote=FALSE);
  }  
}

multiSample<-function(samples, chr){
  ans<-readNG(samples[1],chr=chr);
  for(i in 2:length(samples)){
    x<-readNG(samples[i],chr=chr);
    ans<-cbind(ans, as.vector(x$nimblegen_signal));
  }
  colnames(ans)<-c("position", as.vector(samples));
  ans;
}

# in: data.frame of sample names and file names, and col specs
# out: each file parsed correctly for input to analysis routines and plotting
loadSamples<-function(samples,iPos,iSignal,iChr,ncols,autos=c(1:19),sex=c(TRUE, TRUE), prepend0=FALSE,suffix='.wutxt'){
  for(i in seq(nrow(samples))){
    fname<-samples[i,]$file_name;
    sname<-samples[i,]$sample_name;
    stripNGfile(fname,sname,iPos,iSignal,iChr,ncols,autos,sex,prepend0,suffix);
  }
}
