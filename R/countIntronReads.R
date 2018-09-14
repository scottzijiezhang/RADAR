#' @title countIntronReads
#' @param bamFolder Path to the folder where bam file locates
#' @param binSize The size of consecutive bins to slice the transcripts
#' @param threads The number of threads to use for hyperthreading
#' @param strandToKeep According to library preparation protocol, choose which strand to count. Stranded RNA library usually seq the "ooposite" strand. Small RNA library seq the "same" strand.
#' @param outputDir The directory to save output files
#' @param paired Logical indicating whether the input bam files are from paired end sequencing. Default is FALSE. If using paired end data, the read length will be estimated from the data and only good mate are counted.
#' @export
countIntronReads<-function(
  samplenames,# file name of samples
  gtf, # gtf file used for peak calling
  fragmentLength = 150, 
  bamFolder,
  strandToKeep = "opposite",
  paired = FALSE,
  threads = 1
){

  ## Check bam files
  bamPath.input = paste(bamFolder,"/",samplenames,".input.bam",sep="")
  no.samples = length(samplenames)
  
  if( !all(file.exists(bamPath.input)) ) stop( "input bam file missing!!!" )
  num_bam_files <- length(bamPath.input)
  for (ibam in 1:num_bam_files) {
    inputfile = bamPath.input[ibam]
    if (! file.exists(paste(inputfile,'.bai',sep=""))) {
      print(paste("Stage: index bam file", inputfile))
      indexBam(inputfile)
    }
  }
 
  ## This step removes ambiguous annotations and returns gene model
  cat("Reading gtf file to obtain gene model\nFilter out ambiguous model...\n")
  geneGRList = gtfToGeneModel(gtf) #get the gene model in GRList with only single chromosome and strand.
  cat("Gene model obtained from gtf file...\n")
  
  ## get intronic regions
  txdb <- makeTxDbFromGFF(gtf,format = "gtf")
  Introns <- intronicParts(txdb)
  
  no.genes=length(unique(unlist(Introns$gene_id))) # count number of genes
  
  cat("counting reads for each genes, this step may takes a few hours....\n")
  start_time <- Sys.time()
  registerDoParallel( cores = threads)
  cat(paste("Hyper-thread registered:",getDoParRegistered(),"\n"))
  cat(paste("Using",getDoParWorkers(),"thread(s) to count reads in continuous bins...\n"))
  reads <- foreach(i = 1:no.genes, .combine = rbind, .errorhandling = "remove" ) %dopar%{
    geneName = unique(unlist(Introns$gene_id))[i]
    geneModel = Introns[ unlist(Introns$gene_id) == geneName]
    geneModel <- geneModel[width(geneModel)>20]# filter out too small introns
    
    if(strandToKeep == "opposite"){
      reads.strand <- ifelse(all(strand(geneModel)=="+"),"-","+")# switch strand if necessary
    }else if(strandToKeep == "same"){
      reads.strand <- unique(as.character(strand(geneModel)))
    }else{
      stop("Must provide strand to keep!")
    }
    
    inputCount <- foreach(gr = 1:length(geneModel), .combine = rbind)%do%{
      sapply(bamPath.input, .countBam, Range = geneModel[gr],reads.strand = reads.strand, fragmentLength = fragmentLength, paired = paired  )
    }
    
    rownames(inputCount) <- paste(geneName,1:length(geneModel),sep = ".")
    inputCount
  }
  rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
  end_time <- Sys.time()
  cat(paste("Time used to count reads:",difftime(end_time, start_time, units = "mins"),"mins... \n"))
  colnames(reads) <- samplenames
  
  return(reads)
}


.countBam <-
  function(bam, Range , reads.strand, fragmentLength, paired){
    
    if(paired){  # pair-end sequence
      ba = scanBam(bam, param=ScanBamParam( which=Range, what =c("pos","strand","qwidth","isize"), flag = scanBamFlag(isProperPair = T,isSecondMateRead = F,isFirstMateRea =T) ) )
      ba = data.frame( pos = ba[[1]]$pos, strand = ba[[1]]$strand ,readLength = ba[[1]]$qwidth, isize = ba[[1]]$isize )
      ba = ba[ba$strand == reads.strand, ] ## filter for strand
      ba = ba[abs(ba$isize) <500 , ] # filter for insert size
      ##shift the read pos to the center of the reads
      if(reads.strand == "+"){ba$center = ba$pos + round(abs(ba$isize)/2) }else{ba$center = ba$pos + ba$readLength - round(abs(ba$isize)/2)  }
      ba = ba[ba$center > start(Range), ]
      ba = ba[ba$center < end(Range), ]
      return(nrow(ba))
      
    }else{ # single-end sequence
      ba = scanBam(bam, param=ScanBamParam( which=Range, what =c("pos","strand","qwidth") ) )
      ba = data.frame( pos = ba[[1]]$pos, strand = ba[[1]]$strand, readLength = ba[[1]]$qwidth )
      ba = ba[ba$strand == reads.strand, ] ## filter for strand
      ##shift the read pos to the center of the reads
      if(reads.strand == "+"){ba$pos = ba$center + round(fragmentLength/2) }else{ba$center = ba$pos + ba$readLength - round(fragmentLength/2)  }
      ba = ba[ba$center > start(Range), ]
      ba = ba[ba$center < end(Range), ]
      ba = ba[!is.na(ba$center),]
      return(nrow(ba))
    }
    
  }
