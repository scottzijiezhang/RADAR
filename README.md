# RNADMethyl
Tool sets to analyze MeRIP-seq high throughput data for differential RNA modifications

### Install the R package from Github

Depends: Rsamtools, GenomicFeatures (>= 1.14.5), BH

	install.packages("devtools")
	library(devtools)
	install_github("scottzijiezhang/RNADMethl")
	library("RNADMethl")

### Example to divide gene into 50bp bins and count read
```
RDM <- countReads(
  samplenames, # The program search for paired files for each sample: sample1.input.bam + sample1.m6A.bam
  gtf, # Path to GTF annotation file
  shift = 75, # This is usually half of your RNA fragment size (Usually ~150bp in our lab)
  bamFolder, # Path to the folder where your bam files located
  outputDir=NA, # Path to the folder where you want to save the read count result
  modification = "m6A", # the post fix of IP sample bam file usually sample1.m6A.bam
  binSize = 50, # What's the size of bin you want to slice the transcript 
  strandToKeep = "opposite", # For stranded library protocol, you reads is at the opposite strand of the RNA 
  threads = 1 # number of thread to use 
)
```
### Example to run normalization and filtering of the read count
Normalize the input read count by geometry mean implemented in DESeq2. The ip read count is normalized based on geometry mean of IP efficiency (enrichment) of top-read-count bins.  
```
X <- c(rep("Ctl",8),rep("Case",8))
RDM <- normalizeLibrary( RDM, X )
```
Normalize the IP read count by gene-wise size factor based on INPUT gene level count.
```
RDM <- adjustExprLevel( RDM )
```
Filter out the bins with too little read count (insufficient coverage) and bins not enriched in IP-expriment.  
```
RDM <- filterBins( RDM ,minCountsCutOff = 15)
```

### Example to run PoissonGamma test without covariate
```
RDM <- diffIP(RDM)

## If you want to excude certain sample from the test because you think it is outlier:
RDM <- diffIP(RDM, exclude = "Sample1" )
```
Report the result at FDR 0.1 as bed format.
```
reportPoissonGammaMerge(RDM,cutoff = 0.1)
```

### Example to run PoissonGamma test with covariates
```
X2 <- as.fumeric(c("A","A","A","B","B","A","A","A","A","A","B","A","B","B","A","A"))-1
X3 <- as.fumeric(c("A","A","A","A","A","C","C","C","A","A","A","A","A","A","A","C"))-1
X4 <- as.fumeric(c("M","M","F","M","M","M","F","M","M","M","M","F","M","M","M","M","F"))-1
cov <- cbind(X2,X3,X4)

RDM <- diffIP(RDM,Covariates = cov)
```