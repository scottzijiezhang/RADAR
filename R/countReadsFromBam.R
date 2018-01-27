## a helper function to count reads from bam files
.countReadFromBam <-
  function(bam,which,DNA2RNA,reads.strand,shift,left,sliding){
    ba = scanBam(bam, param=ScanBamParam( which=which, what =c("pos","strand") ) )
    ba = data.frame( pos = ba[[1]]$pos, strand = ba[[1]]$strand )
    ba = ba[ ba$pos > left, ]
    ba = ba[ba$strand == reads.strand, ] ## filter for strand
    ba$pos = DNA2RNA[ba$pos - left] ## convert mapped read pos into RNA position
    ba = ba[ ba$pos > 0, ] ## drop intron reads.
    ##shift the read pos to the center of the reads
    if(reads.strand == "+"){ba$pos = ba$pos + shift}else{ba$pos = ba$pos +50 - shift}
    ##count the reads in the sliding windows
    no.window = length(sliding)
    windowCounts = vector(length = no.window)
    for(j in 1:no.window){
      windowCounts[j]= sum( abs(sliding[j] - ba$pos) <= 25 )
    }
    return(windowCounts)
  }
