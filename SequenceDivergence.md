Using an HMM for detecting regions of sequence divergence between two genomes

wuHMM can utilize sequence divergence to reduce the CNV false positive rate.  To accomplish this, the algorithm needs to be supplied with regions that are deemed 'divergent' from the reference genome.  Below are instructions for one method to estimate such regions given dense genotype calls.

Code:

http://wuhmm.googlecode.com/files/forSeqDiver.R

http://wuhmm.googlecode.com/files/hmm_headSD

http://wuhmm.googlecode.com/files/hmm_footSD

Sample input files:

Genotypes for sample 1: http://wuhmm.googlecode.com/files/samp_7_chr18.txt

Genotypes for sample 2: http://wuhmm.googlecode.com/files/samp_17_chr18.txt


Input genotype files should be tab-delimited text files of the format (NO HEADER):

```
SNP ID          position genotype call    confidence
NES17425611	3008456	A	0.999168
NES17425604	3008870	T	0.998753
NES17425589	3009619	T	0.998575
NES17425591	3010064	C	0.999168
NES17424778	3199629	A	0.999168
NES17424779	3199637	T	0.999093
NES17424780	3200200	A	0.999168
NES17424781	3200567	A	0.999168
NES17424782	3200606	T	0.999168
NES17424783	3201327	C	0.999093
...
```


Launch R from the directory containing forSeqDiver.R, hmm\_headSD, and hmm\_footSD.

```
# Load the functions defined in forSeqDiver.R
> source("forSeqDiver.R")

# now run HMM on chr 18
> assignBlocks(chr=18,slabel=17,"samp_17_chr18.txt","samp_7_chr18.txt", min_states=5)

# by default, output will be written to file 'seqdiver_7_chr18.txt':
#str    stop    num SNPs  state (0=divergent, 1=similar, 2=No calls)
3008456 3488693 273     1
3489295 3517916 27      2
3518260 4594128 1488    1
4594138 4620980 72      0
4622024 4637128 39      1
4637373 4645616 46      0

```