#' @title countReadsCarRNA
#' @description This is the very first function in MeRIP-seq data analysis that initianize a `MeRIP` object. This function takes BAM files of Input/IP library of each samples as input and use given GTF file as gene annotation to divide genes into consecutive bins of user defined size.  
#' @param samplenames The names of each sample (prefix for bam files).
#' @param gtf The gtf format gene annotation file
#' @param fragmentLength The RNA fragment length (insert size of the library).
#' @param modification The modification used to name the BAM files.
#' @param bamFolder Path to the folder where bam file locates
#' @param binSize The size of consecutive bins to slice the transcripts
#' @param nonCodingAnnotation The annotation file or files for non-coding RNA in SAF format. Can be a character or vector of path(s) to the SAF file(s). The default is NA.
#' @param threads The number of threads to use for hyperthreading
#' @param strandToKeep According to library preparation protocol, choose which strand to count. Stranded RNA library usually seq the "ooposite" strand. Small RNA library seq the "same" strand.
#' @param outputDir The directory to save output files
#' @param saveOutput Logical option indicating whether to save output as an RDS file.
#' @param paired Logical indicating whether the input bam files are from paired end sequencing. Default is FALSE. If using paired end data, the read length will be estimated from the data and only good mate are counted.
#' @export
countReadsCarRNA<-function(
  samplenames,# file name of samples
  gtf, # gtf file used for peak calling
  fragmentLength = 150,
  bamFolder,
  nonCodingAnnotation = NA,
  outputDir=NA,
  modification = "m6A",
  binSize = 50,
  strandToKeep = "opposite",
  paired = FALSE,
  threads = 1,
  saveOutput = FALSE
){
  
  ##read bam files
  bamPath.input = paste(bamFolder,"/",samplenames,".input.bam",sep="")
  bamPath.IP = paste(bamFolder,"/",samplenames,".",modification,".bam",sep="")
  no.samples = length(samplenames)
  
  ## Check for missing files and index bam files
  if( !all(file.exists(bamPath.input)) ) stop( "input bam file missing!!!" )
  if( !all(file.exists(bamPath.IP)) ) stop( "IP bam file missing!!!" )
  num_bam_files <- length(bamPath.input)
  for (ibam in 1:num_bam_files) {
    inputfile = bamPath.input[ibam]
    IPfile = bamPath.IP[ibam]
    if (! file.exists(paste(inputfile,'.bai',sep=""))) {
      print(paste("Stage: index bam file", inputfile))
      indexBam(inputfile)
    }
    if (! file.exists(paste(IPfile,'.bai',sep=""))) {
      print(paste("Stage: index bam file", IPfile))
      indexBam(IPfile)
    }
  }
  
  
  ## This step removes ambiguous annotations and returns gene model
  cat("Reading gtf file to obtain gene model\nFilter out ambiguous model...\n")
  geneGRListFromGTF = gtfToGeneModel(gtf) #get the gene model in GRList with only single chromosome and strand.
  cat("Gene model obtained from gtf file...\n")
  geneGRList.wholeGene <- range(geneGRListFromGTF) # collpase exons into one transcript
  
  ## Build model for non-coding RNAs if annotation is supplied.
  if( any( ! is.na( nonCodingAnnotation) ) ){
    cat( "Non-coding RNA model obtained from SAF file...\n" )
    if( length(nonCodingAnnotation) ==1 ){
      nonCodingAnnotation.saf <- read.table(nonCodingAnnotation, sep = "\t", header = TRUE, col.names = c("GeneID", "Chr", "Start", "End","Strand") )
    }else{
      nonCodingAnnotation.saf <- foreach( ann = nonCodingAnnotation, .combine = rbind )%do%{ read.table( ann, sep = "\t", header = TRUE, col.names = c("GeneID", "Chr", "Start", "End","Strand") ) }
    }
    
    geneGRList.nonCoding <- makeGRangesListFromDataFrame(nonCodingAnnotation.saf, names.field = "GeneID")
    names(geneGRList.nonCoding) <- nonCodingAnnotation.saf$GeneID
    
    geneGRListCombine <- c( geneGRList.wholeGene, geneGRList.nonCoding )
  }else{
    geneGRListCombine <- geneGRList.wholeGene
  }
  
  ## Check BAM headers and remove chr in geneModel that is not in BAM file. 
  bamHeader <- scanBamHeader(bamPath.input, what=c("targets") )
  seqLevels <- unique( unlist( lapply( bamHeader, function(x) names( x$targets) ) ) )
  geneGRListCombine <- geneGRListCombine[ unlist( runValue( seqnames( geneGRListCombine ) ) ) %in% seqLevels ]
  
  
  ## Process the gene region first.
  ## Count reads in transcripts. Because we are interested in carRNA, we count reads on pre-mRNA instead of mRNA. 
  no.genes=length(geneGRListCombine)## define number of genes
  
  ## Divide annotations into batches for counting reads
  batchNum <-  max( round( no.genes/2e4 ), 1 ) # start with at least 1 batch
  batchPoints <- round( seq( 1, no.genes, length.out = batchNum+1  ) ) # Stop points of each batch
  batches <- foreach( j = 1:(length(batchPoints) -1 ) )%do%{ return(batchPoints[j]:batchPoints[j+1] ) }
  
  cat(paste("Dividing annotations into", batchNum,"batches for read count quantification... \n"))
  
  
  for( iBatch in 1:batchNum ){
    
    cat(paste0( "Counting reads for each gene in batch # ",iBatch, ", this step may takes a few hours...\n") )
    
    tmpReads <- .CountReadsBatch(geneGRList = geneGRListCombine[ batches[[iBatch]] ],
                                 binSize = binSize,
                                 bamPath.IP = bamPath.IP,
                                 bamPath.input = bamPath.input,
                                 strandToKeep = strandToKeep,
                                 fragmentLength = fragmentLength,
                                 paired = paired,
                                 threads = threads
    )
    
    cat(paste0( "Finished quantification for batch # ",iBatch,"; ",batchNum - iBatch, " batches to be processed...\n\n") )
    
    ## filter out zero count genes
    rowTotoal <- rowSums(tmpReads) > 0
    tmpGeneName <- unlist( lapply( strsplit( rownames(tmpReads) , ","), function(x) x[1] ) )
    tmpZeroGene <-  tapply(rowTotoal , tmpGeneName, sum) 
    if( any(tmpZeroGene == 0) ){
      tmpReads <- tmpReads[ tmpGeneName %in% names( tmpZeroGene[tmpZeroGene >= 0] ), ]
      
    }
    
    ## save the current data
    eval(parse(text = paste0("Reads_",iBatch," <- tmpReads ") ) )
    
  }
  
  cat(paste0( "Finished reads quantification. Combining read counts from batch runs...\n") )
  
  ## combine all reads
  eval(parse(text = paste0("reads <- rbind(",paste0( "Reads_",1:batchNum,collapse = ",") , ") ") ) )
  
  ## remove temp read count tables
  eval(parse(text = paste0("rm(",paste0( "Reads_",1:batchNum,collapse = ",") , ") ") ) )
  
  colnames(reads) <- c(paste(samplenames,"input",sep = "-"),paste(samplenames,"IP",sep = "-"))
  
  
  data.out <- MeRIP(reads = reads, binSize = binSize, gtf = gtf, geneModel = geneGRListCombine, bamPath.input = bamPath.input, bamPath.ip = bamPath.IP, samplenames = samplenames, mode = "carRNA")
  if(saveOutput){
    ## create output directory
    dir.create(outputDir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(data.out,paste0(outputDir,"/MeRIP_readCounts.RDS"))
  }
  
  
  return(data.out)
}


## A helper function to process a batch of annotation features.
.CountReadsBatch <- function(geneGRList,
                             binSize,
                             bamPath.IP,
                             bamPath.input,
                             strandToKeep,
                             fragmentLength,
                             paired,
                             threads
                             ){
  
  
  start_time <- Sys.time()
  registerDoParallel( cores = threads)
  cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
  cat(paste("Using",getDoParWorkers(),"thread(s) to count reads in continuous bins...\n"))
  reads <- foreach(i = 1:length(geneGRList), .combine = rbind) %dopar%{
    
    geneName = names(geneGRList)[i]
    geneModel =geneGRList[geneName][[1]] 
    
    # DNA location to gene location conversion
    df.geneModel= as.data.frame(geneModel) ##data frame of gene model
    dna.range = as.data.frame(range(geneModel))
    df.geneModel$end = df.geneModel$end - dna.range$start + 1
    df.geneModel$start = df.geneModel$start - dna.range$start + 1
    DNA2RNA = 1:dna.range$width
    #DNA2RNA = rep(0,dna.range$end - dna.range$start +1)
    #no.exon = dim(df.geneModel)[1]
    #for (j in 1:no.exon){DNA2RNA[df.geneModel$start[j]:df.geneModel$end[j]]=1}
    #exon.length = sum(DNA2RNA)
    #DNA2RNA=cumsum(DNA2RNA)*DNA2RNA
    geneLength = dna.range$width
    
    ## skip any gene with smaller than 200bp transcript
    if(geneLength < 200) {return(NULL)}
    
    
    ## switch strand because stranded RNA library protocol sequence reverse strand
    if(strandToKeep == "opposite"){
      reads.strand = character()
      if(dna.range$strand == "+"){reads.strand = "-"}else if(dna.range$strand == "-"){reads.strand = "+"}else{reads.strand = "*"} ## switch strand on RNA reads for Truseq protocol
    }else if(strandToKeep == "same"){
      reads.strand = as.character(dna.range$strand)
    }else{
      cat("Currently RADAR only support strand specific RNA-seq data.\nCounting reads at opposite strand by defalt...\n")
      reads.strand = character()
      if(dna.range$strand == "+"){reads.strand = "-"}else if(dna.range$strand == "-"){reads.strand = "+"}else{reads.strand = "*"}
    }
    
    #create start points of continuous window
    if(exon.length <= binSize){
      slidingStart = 1
    }else{
      ## use the 3' end terminal bin as a elastic-size bin
      if(dna.range$strand == "+"){
        slidingStart = seq(from = 1, to = ( exon.length - binSize - exon.length %% binSize + 1) , length.out = floor(exon.length/binSize) ) 
      }else{ # make the first bin elastic bin if a gene is on reverse strand
        slidingStart = c(1, seq(from = binSize + exon.length %% binSize + 1, to = ( exon.length - binSize + 1) , length.out = floor(exon.length/binSize) - 1 )  )
      }
    }
    
    #mapping$chr = as.character(dna.range$seqnames)
    #mapping$strand = as.character(dna.range$strand)
    #rownames(mapping) = paste(geneName,slidingStart,sep = ",")
    #geneRNA2DNA= rbind(geneRNA2DNA,mapping[c("chr","start","end","strand")])
    
    #count reads in all samples
    ba.IP = sapply(bamPath.IP,.countReadFromBam,which = range(geneModel),reads.strand = reads.strand,DNA2RNA = DNA2RNA,fragmentLength=fragmentLength,left=dna.range$start,sliding = slidingStart, binSize = binSize, paired = paired)
    ba.input = sapply(bamPath.input,.countReadFromBam,which = range(geneModel),reads.strand = reads.strand,DNA2RNA = DNA2RNA,fragmentLength=fragmentLength,left=dna.range$start,sliding = slidingStart, binSize = binSize, paired = paired)
    
    if(is.vector(ba.IP) ){# if there is only one window for this gene, make it a matrix to avoid bug
      ba.IP = matrix(ba.IP, nrow = 1)
      ba.input = matrix( ba.input, nrow = 1 )
    }
    ba.counts <- cbind(ba.input,ba.IP)
    rownames(ba.counts) <-  paste(geneName,slidingStart,sep = ",")
    
    ba.counts
  }
  rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
  end_time <- Sys.time()
  cat(paste("Time used to count reads in this batch:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
  
  return(reads)
  
}










