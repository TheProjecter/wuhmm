First, download and unzip [wuHMM](http://wuhmm.googlecode.com/files/wuHMM1_1.zip) or follow [these instructions](http://code.google.com/p/wuhmm/source/checkout) for obtaining the latest wuHMM source code.  The example R session below describes the process of loading data files, using wuHMM, and plotting the results.  The data files used in the examples below are available for download from [here](http://compbio.s3.amazonaws.com/demoData.zip).  Download and decompress (unzip) the compressed demo data file in the same directory that you placed wuHMM.  This data is also available at GEO under accession [GSE5805](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5805).  Download the sample text file that contains the demo sample data descriptions [here](http://wuhmm.googlecode.com/files/project.txt).

Launch R.   The following R commands will get you started reading CGH data, detecting, genotyping, and plotting CNV predictions.  Lines beginning with '>' are commands, '#' are comments.

### Loading data ###

```
# Read sample file containing sample and file names.
> samples<-read.table("project.txt", sep="\t",quote='', header=TRUE, as.is=TRUE)

# Prepare demo data files for analysis.
# The demo data files are from NimbleGen, however other data types should be understandable by the algorithm, provided that the file format is specified correctly.

# You need to specify the column that contains the following data:
# Probe position: iPos = 5
# Chromosome: iChr = 3
# Log2(ratio), intensity, read count: iSignal = 14

# Other parameters:
# Total number of columns in file: cols = 14
# Chromosome names:
#   numeric chrs: autos
#   sex chrs: boolean vector
# Chromosome format.  Double-digit (i.e. Chr01: prepend0=TRUe), or single (Chr1)?  
> system.time(loadSamples(samples,iPos=5,iSignal=14,iChr=3,ncols=14,autos=c(1:19),sex=c(TRUE,FALSE), prepend0=0))
[1] 50.14  3.31 56.50  0.00  0.00

```

### Running wuHMM ###
```
# Decode (or segment) all samples, all chromosomes (except ChrY).  wuHMM will output the sample name as it is decoding.
# wuHMM parameters, default parameters listed:
#    sample_labels (required): must be same as samples$sample_name used in loadSamples.
#    chrs: Chromosomes to segment
#    nstay_seed=4:  Number of probes required to 'seed' a region prior to HMM processing.    Increasing this number will reduce the FPR AND sensitivity.
# Other parameters can be found in funs.R.
> system.time(cnvsAll<-decodeAndPerm(samples$sample_name,chrs=c(1:19,'x'), nstay_seed=5))
AJ
Balb_cByJ
C57BL_6J
C58_J
DBA_2J
[1]  65.07  72.22 427.89 162.55 118.38

# How many CNV predictions?
> dim(cnvsAll)
[1] 148   9

# Filter out low-scoring calls
# Parameters: cnvs, score thresholds for (1) losses, and (2) gains
> cnvsGood<-filterCalls(cnvsAll, 1.5,1)
> dim(cnvsGood)
[1] 87  9

# List first 3 CNVs order by str
> cnvsGood[order(cnvsGood$str),][1:3,]
        str     stp num_probes   mean_sig    score     noise median_sig chr name
132 6110655 6165058          9  0.5220000 1.021709 0.1995345      0.465  17 Balb_cByJ
28  6413122 8404697         47 -0.5714043 2.456394 0.2677412     -0.638   7 Balb_cByJ
142 9109176 9545727         40 -0.7248250 2.596971 0.2780237     -0.704   7 C58_J

# Plot a CNV.
# Parameters (required): sample_name, cnv
# Other parameters, margin size, y-axis boundaries, etc, can be found in plotting.R.
# plotReg will display all cnvs supplied if they are within margin distance of the first CNV.

> plotReg('Balb_cByJ',cnvsGood['133',], marg=50)
```
![http://wuhmm.googlecode.com/files/cnvMarg50.gif](http://wuhmm.googlecode.com/files/cnvMarg50.gif)
```
# Now zoom in (default setting):
> plotReg('Balb_cByJ',cnvsGood['133',])
```
![http://wuhmm.googlecode.com/files/cnv.gif](http://wuhmm.googlecode.com/files/cnv.gif)


### Merge and Genotype CNVRs ###
```
# Merge CNVs into CNV-regions (CNVRs).  Assign each CNVR an identifier starting with strID.
# makeCNVRs returns a list of the following elements:
#  (1) CNVRs
#  (2) List of CNV tables.
> cnvrs<-w_makeCNVRs(cnvsGood, 1)
> dim(cnvrs[[1]])
[1] 58  6

# list 3 CNVRs
> cnvrs[[1]][1:3,]
  cnvrid chr       str       stp length ave_score
1      1   1 100982630 101050290  67661     2.942
2      2   1 171447273 171508979  61707     2.599
3      3   5  14292154  14539614 247461     2.052

# Assign an inferred genotype to each sample based on the optimal partitioning of the data .
# Paramters: sample_names, list of CNV tables and CNVRs from running w_makeCNVRs.
> cnvrGenos<-w_genotype(samples$sample_name, cnvrs[[2]], cnvrs[[1]])

# Plot all samples together, highlighting samples with 'abnormal' genotypes.
# First a CNVR comprised of gains
> w_plotCNVR(cnvrGenos, 3)
```
![http://wuhmm.googlecode.com/files/cnvr3.gif](http://wuhmm.googlecode.com/files/cnvr3.gif)
```
# Now a CNVR comprised of losses
> w_plotCNVR(cnvrGenos, 26)
```
![http://wuhmm.googlecode.com/files/cnvr26.gif](http://wuhmm.googlecode.com/files/cnvr26.gif)