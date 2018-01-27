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
countReads(
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
