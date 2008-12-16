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

decodeAndPerm<-function(sample_labels,
                        nperms=0,
                        nclusts=7,
                        diver=FALSE,          
                        chrs=c(1:19),
                        nstates=3,
                        nstay_seed=4,
                        scale=5,
                        fact=0,
                        files=c("hmm_head","hmm_body","hmm_foot"),
                        em=1,
                        iterations=100,
                        lthresh=.1,
                        marg=5,
                        nstay_pures=2){

  allCalls<-data.frame();
  nCNVs<-0;
  for(i in seq(length(sample_labels))){
    sample_label<-sample_labels[i];
    for(chr in chrs){
      fname_cnvs<-paste("cnvs_",sample_label,".csv",sep='');
      dfX<-readNG(sample_label,chr=chr);

      # Cluster log2-ratios 
      chr2<-paste("chr",chr,sep='');
      rownames(dfX)<-rank(dfX$position, ties.method="rand");
      ng<-addClustRank(dfX, nclusts, diver);
      calls<-data.frame();
      if(nrow(ng)>0){
        # decode unshuffled data
        calls<-decode(ng=ng,
                      samp=sample_label,
                      nclust=nclusts,
                      norm=norm,
                      abnorm=abnorm,
                      joinNorm=joinNorm,
                      joinAbnorm=joinAbnorm,
                      nstay_pures=nstay_pures,
                      iterations=iterations,
                      lthresh=lthresh,
                      nstay_seed=nstay_seed,
                      marg=marg,
                      scale=scale,
                      fact=fact);
        # decode shuffled data
        if(nperms>0 & nrow(calls)>0 ){
          rscores<-vector(length=nperms);
          for(perm_n in seq(nperms)){
            ng2<-ng;
            rows<-sample(seq(nrow(ng)),nrow(ng));
            ng2$position<-ng[rows,]$position;
            ng2<-ng2[order(ng2$position),];
            rCalls<-decode(ng2,
                           samp=samp,
                           nclust=nclusts,
                           norm=norm,
                           abnorm=abnorm,
                           joinNorm=joinNorm,
                           joinAbnorm=joinAbnorm,
                           nstay_pures=nstay_pures,
                           iterations=iterations,
                           lthresh=lthresh,
                           nstay_seed=nstay_seed,
                           marg=marg,
                           scale=scale,
                           fact=fact);
            if(nrow(rCalls)>0){
              rscores[perm_n]<-max(rCalls$score);
            }
            else{
              rscores[perm_n]<-0;
            }
          }
          # compute p-values
          pvals<-getpvals(calls$score,rscores);
          calls<-cbind(calls,pvals);
          colnames(calls)[ncol(calls)]<-"pval";
        }        
      }
      calls<-cbind(calls,data.frame(chr=rep(chr, nrow(calls))));
      calls<-cbind(calls,data.frame(name=rep(sample_label, nrow(calls))));
      allCalls<-rbind(allCalls, calls);
    }
  }
  if(nperms>0){
    colnames(allCalls)[ncol(allCalls)]<-"pval";
  }
  allCalls;
}
  

decode<-function(ng,
                 samp,
                 nclust=7,
                 nstates=3,
                 nstay_seed=3,
                 scale=5,
                 fact=1,
                 files=c("hmm_head","hmm_body","hmm_foot"),
                 em=1,
                 iterations=100,
                 lthresh=1,
                 marg=5,
                 nstay_pures=3:10,
                 norm=0.6,
                 abnorm=0.93,
                 joinNorm=0.75,
                 joinAbnorm=0.9){

  ab_states<-1:(nstates-1);
  cnvsAll<-data.frame();

  eseeds<-seedPython(ng, nclust, nstay_seed);

  if(nrow(eseeds)>0){
    eseeds$num_probes<-countSeedProbes(ng, eseeds);
    seeds<-mergeSeeds(expandSeeds(ng,eseeds,marg));

    # 'finish the ends'
    nseeds<-nrow(seeds);
    seeds$num_probes<-countSeedProbes(ng, seeds);
    seeds<-mergeSeeds(expandSeeds(ng,seeds,0.25));

    if(nrow(seeds)!=nseeds){
      seeds$num_probes<-countSeedProbes(ng, seeds);
      seeds<-mergeSeeds(expandSeeds(ng,seeds,0.25));
    }
    
    cat(samp,"\n");
    for(i in seq(nrow(seeds))){
      seed<-seeds[i,];
      #cat("seed ",i, seed$str,"-",seed$stp,"\n");
      ng2<-ng[ng$position>=seed$str & ng$position<=seed$stp,];
      phits<-decodeOpt(ng=ng2,
                       samp=samp,
                       nclusts=nclust,
                       norm,
                       abnorm=abnorm,
                       joinNorm,
                       joinAbnorm,
                       nstay_pures=nstay_pures,
                       iterations=iterations,
                       lthresh=lthresh);
      # test for overlap with eseeds
      if(nrow(phits)>0){
        for(j in seq(nrow(phits))){
          phit<-phits[j,];
          for(z in seq(nrow(eseeds))){
            if(overlap(eseeds[z,],phit)){
              cnvsAll<-rbind(cnvsAll, phit);
              break;
            }
          }
        }
      }
    }
    if(nrow(cnvsAll)>0){
      scs<-scores3(cnvs=cnvsAll,scale=scale,ng=ng, fact=fact);
      cnvsAll<-cbind(cnvsAll, scs[,1]);
      colnames(cnvsAll)[ncol(cnvsAll)]<-"score";
      cnvsAll<-cbind(cnvsAll, scs[,2]);
      colnames(cnvsAll)[ncol(cnvsAll)]<-"noise";
      cnvsAll<-cbind(cnvsAll, scs[,3]);
      colnames(cnvsAll)[ncol(cnvsAll)]<-"median_sig";

    }
  }
  
  cnvsAll;
}

seedPython<-function(ng, nclusts=9, min_seed=3){
  write.table(ng, file="ngForSeeds.csv", sep="\t", col.names=FALSE, row.names=FALSE,quote=FALSE);
  cmd<-paste("python seedPython.py"," ",nclusts," ",min_seed,sep='');
  #cat(cmd,"\n");
  system(cmd, intern=TRUE);
  seeds<-read.table(file="mySeeds.csv",sep='\t', header=TRUE);
  seeds;
}


decodeOpt<-function(ng,
                    samp,
                    nclusts,
                    norm,
                    abnorm,
                    joinNorm,
                    joinAbnorm,
                    nstates=3,
                    nstay_pures=3:10,
                    files=c("hmm_head","hmm_body","hmm_foot"),
                    em=1,
                    iterations=100,
                    lthresh=1){

  nstay_pure<-nstay_pures[1];
  lt_ts<-data.frame(ab_self=c(0.704,0.554, 0.304),
                    ab_norm=c(0.244, 0.394, 0.644),
                     ab_joiner=c(0.001, 0.001, 0.001));
#                    ab_joiner=c(0.05, 0.05, 0.05));
#  jts<-data.frame(j_self=c(0.001), j_same=c(0.999));

  jts<-data.frame(j_self=c(0.1, 0.25), j_same=c(0.9, 0.75));
  abnorm<-c(0.999);
  norms<-c(0.4, 0.7);
  joinAbnorm <- c(0.4);
  joinNorm <- c(0.7);
  
  min_ll<-1e8;
  best_nstay<-nstay_pures[1];

  best_lt<-1;
  best_jt<-1;
  best_norm<-1;

  fname<-paste("clust_",samp,".csv",sep='');
  write.table(ng, file=fname, sep="\t", col.names=FALSE, row.names=FALSE,quote=FALSE);
  
  for(i in seq(nrow(lt_ts))){
    for(j in seq(nrow(jts))){
      for(a in seq(length(norms))){
#  for(nstay_pure in nstay_pures){
        ab_states<-1:(nstates-1);
        fname_cnvs<-paste("cnvs_",samp,".csv",sep='');
        cnvsAll<-data.frame();
        sink(files[2]);
        write_emis(nclusts,nstay_pure,norm=norms[a], abnorm=abnorm, joinNorm, joinAbnorm);
        write_trans(nstay_pure,lt=lt_ts[i,],jt=jts[j,],seed=FALSE);
        write_dict(nstay_pure);
        write_pi(nstay_pure);
        cat("\nlthresh=",lthresh,";\n",sep='');
        sink();
        hmm_body(nclusts, iterations=iterations,files[2],samp,fname_cnvs,ab_states=ab_states);
        pyth_script<-hmm_build(samp, c(files[1:2],"hmm_foot_sel"));
        
        cmd<-paste("python ",pyth_script,sep='');
        ll<-as.numeric(system(cmd, intern=TRUE));
        if( (-1 * ll) < min_ll){
          min_ll<- -1 * ll;
          best_lt<-i;
          best_jt<-j;
          best_norm<-a;
        }

        #cleanup
        #cmd<-paste("rm ",fname_cnvs,sep='');
        #system(cmd);        
        cmd<-paste("rm ", pyth_script,sep='');
        system(cmd);
        
      }
    }
  }
  # cat("best self:",lt_ts[best_lt,]$ab_self,"\n");
  # cat(norms[best_norm],"\n");
  # cat(jts[best_jt,]$j_self,"\n");
  nstay_pure<-best_nstay;
  ab_states<-1:(nstates-1);
  fname_cnvs<-paste("cnvs_",samp,".csv",sep='');
  cnvsAll<-data.frame();
  sink(files[2]);
  write_emis(nclusts,nstay_pure,norm=norms[best_norm], abnorm=abnorm, joinNorm, joinAbnorm);
  write_trans(nstay_pure,lt=lt_ts[best_lt,],jt=jts[best_jt,],seed=FALSE);
  write_dict(nstay_pure);
  write_pi(nstay_pure);
  cat("\nlthresh=",lthresh,";\n",sep='');
  sink();
  hmm_body(nclusts, iterations=iterations,files[2],samp,fname_cnvs,ab_states=ab_states);
  pyth_script<-hmm_build(samp, files);
  
  cmd<-paste("python ",pyth_script,sep='');
  system(cmd,intern=TRUE);
  cnvs<-read.table(file=fname_cnvs,sep='\t', header=TRUE);

  #cleanup
  cmd<-paste("rm ",fname_cnvs,sep='');
  system(cmd);
  cmd<-paste("rm ",fname,sep='');
  system(cmd);
  cmd<-paste("rm ", pyth_script,sep='');
  system(cmd);
  
  cnvs;
}

# returns name of python script created
# files: 1=header, 2=model parameters, 3=footer
hmm_build<-function(sample_label,files=c("hmm_head","hmm_body","hmm_foot")){
  #cat(files,"\n")
  fname<-paste("hmm_",sample_label,".py", sep='');
  cmd<-paste("cat ",files[1]," ",files[2]," ",files[3] ," > ",fname,sep='');
  system(cmd);
  fname;
}

hmm_body<-function(nclusts,
                   iterations,
                   filen,
                   sample_label,
                   fname_cnvs,
                   ab_states=1:2){
  sink(filen,append=TRUE);
  cat("\niterations=",iterations,";\n",sep='');
  cat("sigma = IntegerRange(1,",nclusts+1,")\n",sep='')
  cat("model = ghmm.HMMFromMatrices(sigma,ghmm.DiscreteDistribution(sigma), A, B, pi);\n");
  cat("i_alph = sigma;\n");
  cat("clusts = list();\npositions = list();\nsignals = list();\n");
  cat("\nsample_label='",sample_label,"';\n",sep='');
  cat("\nfname_cnv=\"",fname_cnvs,"\"\n",sep='');
  cat("clustsfile = open(\"./clust_",sample_label,".csv\",'r')\n",sep='');
  cat("ab_states = [");  
  for(i in ab_states[1:length(ab_states)-1]){
    cat(i,",",sep='');
  }
  cat(ab_states[length(ab_states)],"];\n",sep='');
  sink();
}

decodeGenome<-function(dfs,
                       sample_labels,
                       binSize=20,
                       nperms=0,
                       nclusts=7,
                       diver=FALSE,          
                       chrs=c(1:19),
                       nstates=3,
                       nstay_seed=2,
                       scale=5,
                       fact=0,
                       files=c("hmm_head","hmm_body","hmm_foot"),
                       em=1,
                       iterations=100,
                       lthresh=.1,
                       marg=5,
                       nstay_pures=2){


  allCalls<-data.frame();
  nCNVs<-0;
  for(i in seq(length(dfs))){
    df<-dfs[[i]];
    sample_label<-sample_labels[i];

    # get average values
    fname_cnvs<-paste("genome_cnvs_",sample_label,".csv",sep='');

    df<-df[order(df$chr, df$position),];
    rownames(df)<-1:nrow(df);
    df<-averageGenome(df, binSize);
    ng<-addClustRank(df, nclusts, diver);
    calls<-data.frame();
    if(nrow(ng)>0){
        # decode unshuffled data
      calls<-decode(ng=ng,
                    samp=sample_label,
                    nclust=nclusts,
                    norm=norm,
                    abnorm=abnorm,
                    joinNorm=joinNorm,
                    joinAbnorm=joinAbnorm,
                    nstay_pures=nstay_pures,
                    iterations=iterations,
                    lthresh=lthresh,
                    nstay_seed=nstay_seed,
                    marg=marg,
                    scale=scale,
                    fact=fact);
        # decode shuffled data
      if(nperms>0 & nrow(calls)>0 ){
        rscores<-vector(length=nperms);
        for(perm_n in seq(nperms)){
          ng2<-ng;
          rows<-sample(seq(nrow(ng)),nrow(ng));
          ng2$position<-ng[rows,]$position;
          ng2<-ng2[order(ng2$position),];
          rCalls<-decode(ng2,
                         samp=samp,
                         nclust=nclusts,
                         norm=norm,
                         abnorm=abnorm,
                         joinNorm=joinNorm,
                         joinAbnorm=joinAbnorm,
                         nstay_pures=nstay_pures,
                         iterations=iterations,
                         lthresh=lthresh,
                         nstay_seed=nstay_seed,
                         marg=marg,
                         scale=scale,
                         fact=fact);
          if(nrow(rCalls)>0){
            rscores[perm_n]<-max(rCalls$score);
          }
          else{
            rscores[perm_n]<-0;
          }
        }
          # compute p-values
        pvals<-getpvals(calls$score,rscores);
        calls<-cbind(calls,pvals);
        colnames(calls)[ncol(calls)]<-"pval";
      }        
    }
    allCalls<-rbind(allCalls, calls);
  }
  if(nperms>0){
    colnames(allCalls)[ncol(allCalls)]<-"pval";
  }
  allCalls;
}
  
reduceNG<-function(df, nprobes){
  chrs<-unique(df$seq_id);
  dfs<-list();

  for(chr in chrs){
    temp<-data.frame();
    ng<-df[df$seq_id==chr,];
    sample_label<-ng[1,]$sample_label;
    
    npoints<-floor(nrow(ng)/nprobes);
    
    vals<-vector(length=npoints);    
    position<-vector(length=npoints);    
    probe_id<-vector(length=npoints);
    
    beg<-1
    end<-nprobes;
    for(i in 1:npoints){
      x<-ng[beg:end,];
      vals[i]<-mean(x$nimblegen_signal);
      beg<-end+1;
      if(i==(npoints-1)){
        end<-nrow(ng);
      }
      else{
        end<-end+nprobes;
      }
      position[i]<-median(x$position);
      probe_id[i]<-x$probe_id[floor( ((end-beg)/2) + beg )];
    }  
    temp<-data.frame(sample_label=rep(sample_label, nrow(ng)),
                     seq_id=rep(chr, nrow(ng)),
                     probe_id=probe_id,
                     position=position,
                     nimblegen_signal=nimblegen_signal,
                     state=rep(0,nrow(ng)));
    dfs<-append(dfs, temp);
  }
  dfs;
}
  
