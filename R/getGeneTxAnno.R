#' @title getGeneTxAnno
getGeneTxAnno <- function(geneGRList,
                          geneNames
                          ){

  geneName = names(geneGRList)[i]
  geneModel =reduce( geneGRList[geneName][[1]] )## merge overlapping exons

  # DNA location to gene location conversion
  df.geneModel= as.data.frame(geneModel) ##data frame of gene model
  dna.range = as.data.frame(range(geneModel))
  df.geneModel$end = df.geneModel$end - dna.range$start + 1
  df.geneModel$start = df.geneModel$start - dna.range$start + 1
  DNA2RNA = rep(0,dna.range$end - dna.range$start +1)
  no.exon = dim(df.geneModel)[1]
  for (j in 1:no.exon){DNA2RNA[df.geneModel$start[j]:df.geneModel$end[j]]=1}
  exon.length = sum(DNA2RNA)
  DNA2RNA=cumsum(DNA2RNA)*DNA2RNA

  ## skip any gene with smaller than 200bp transcript
  if(exon.length < 200) {next}

  #creat a corresponding map from RNA to DNA
  RNA2DNA = 1:exon.length
  pointer = 1
  for (j in 1:no.exon){
    RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1) ]= RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1)] + df.geneModel$start[j] -pointer
    pointer = pointer + df.geneModel$width[j]
  }
  RNA2DNA = RNA2DNA + dna.range$start -1 #back to chromosome coordinates

  #creat center points of continuous window
  if(exon.length <= binSize){
    slidingStart= exon.length/2
    mapping = data.frame(start = RNA2DNA[slidingStart-exon.length/2+1], end = RNA2DNA[slidingStart + exon.length/2]  )
  }else{
    slidingStart= round(seq(from = binSize/2, to = (exon.length - binSize/2), length.out = ceiling(exon.length/binSize) ) )
    mapping = data.frame(start = RNA2DNA[slidingStart - binSize/2 +1], end = RNA2DNA[slidingStart + binSize/2 ]  )
  }

  mapping$chr = as.character(dna.range$seqnames)
  mapping$strand = as.character(dna.range$strand)
  rownames(mapping) = paste(geneName,slidingStart,sep = ",")
  geneRNA2DNA= rbind(geneRNA2DNA,mapping[c("chr","start","end","strand")])


}


.getGeneBins <- function(geneGRList,geneNames,binSize){

  registerDoParallel(cores = 6)
  geneBins <-foreach(i = 1:length(geneNames), .combine = rbind)%dopar%{
    geneName = geneNames[i]
    geneModel =reduce( geneGRList[geneName][[1]] )## merge overlapping exons

    # DNA location to gene location conversion
    df.geneModel= as.data.frame(geneModel) ##data frame of gene model
    dna.range = as.data.frame(range(geneModel))
    df.geneModel$end = df.geneModel$end - dna.range$start + 1
    df.geneModel$start = df.geneModel$start - dna.range$start + 1
    DNA2RNA = rep(0,dna.range$end - dna.range$start +1)
    no.exon = dim(df.geneModel)[1]
    for (j in 1:no.exon){DNA2RNA[df.geneModel$start[j]:df.geneModel$end[j]]=1}
    exon.length = sum(DNA2RNA)
    #creat a corresponding map from RNA to DNA
    RNA2DNA = 1:exon.length
    pointer = 1
    for (j in 1:no.exon){
      RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1) ]= RNA2DNA[pointer:(pointer+df.geneModel$width[j]-1)] + df.geneModel$start[j] -pointer
      pointer = pointer + df.geneModel$width[j]
    }
    RNA2DNA = RNA2DNA + dna.range$start -1 #back to chromosome coordinates
    #creat center points of continuous window
    if(exon.length <= binSize){
      slidingStart= exon.length/2
      mapping = data.frame(start = RNA2DNA[slidingStart-exon.length/2+1], end = RNA2DNA[slidingStart + exon.length/2]  )
      }else{
      slidingStart= round(seq(from = binSize/2, to = (exon.length - binSize/2), length.out = ceiling(exon.length/binSize) ) )
      mapping = data.frame(start = RNA2DNA[slidingStart - binSize/2 +1], end = RNA2DNA[slidingStart + binSize/2 ]  )
      }
    mapping$chr = as.character(dna.range$seqnames)
    mapping$strand = as.character(dna.range$strand)
    cbind(data.frame(geneName,slidingStart),mapping)
  }
  rm(list=ls(name=foreach:::.foreachGlobals), pos=foreach:::.foreachGlobals)
  return(geneBins)
}


