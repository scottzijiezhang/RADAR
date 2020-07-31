## a helper function to count reads from bam files
.countReadFromBam <-
  function(bam,which,DNA2RNA,reads.strand,fragmentLength,left,sliding, binSize ,paired){
    
    if(paired){  # pair-end sequence
      ba = scanBam(bam, param=ScanBamParam( which=which, what =c("pos","strand","qwidth","isize"), flag = scanBamFlag(isProperPair = T,isSecondMateRead = F) ) )
      ba = data.frame( pos = ba[[1]]$pos, strand = ba[[1]]$strand ,readLength = ba[[1]]$qwidth, isize = ba[[1]]$isize )
      ba = ba[ ba$pos > left, ]
      ba = ba[ba$strand == reads.strand | reads.strand == "*", ] ## filter for strand
      ba$pos = DNA2RNA[ba$pos - left] ## convert mapped read pos into RNA position
      ba = ba[ ba$pos > 0, ] ## drop intron reads.
      ##shift the read pos to the center of the reads
      if(reads.strand == "+"){ba$pos = ba$pos + round(abs(ba$isize)/2) }else{ba$pos = ba$pos + ba$readLength - round(abs(ba$isize)/2)  }
      ##count the reads in the sliding windows
      no.window = length(sliding)
      windowCounts = vector(length = no.window)
      for(j in 1:no.window){
        windowCounts[j]= sum( abs(sliding[j] - ba$pos) <= binSize/2 )
      }
      
    }else{ # single-end sequence
      ba = scanBam(bam, param=ScanBamParam( which=which, what =c("pos","strand","qwidth") ) )
      ba = data.frame( pos = ba[[1]]$pos, strand = ba[[1]]$strand, readLength = ba[[1]]$qwidth )
      ba = ba[ ba$pos > left, ]
      ba = ba[ba$strand == reads.strand | reads.strand == "*", ] ## filter for strand
      ba$pos = DNA2RNA[ba$pos - left] ## convert mapped read pos into RNA position
      ba = ba[ ba$pos > 0, ] ## drop intron reads.
      ##shift the read pos to the center of the reads
      if(reads.strand == "+"){ba$pos = ba$pos + round(fragmentLength/2) }else{ba$pos = ba$pos + ba$readLength - round(fragmentLength/2)  }
      ##count the reads in the sliding windows
      no.window = length(sliding)
      windowCounts = vector(length = no.window)
      for(j in 1:no.window){
        windowCounts[j]= sum( abs(sliding[j] - ba$pos) <= binSize/2 )
      }
    }
    
    return(windowCounts)
  }
