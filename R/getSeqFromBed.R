getSeqFromBed <- function(inputFilePath, header=FALSE, BSgenomeName, upstream = 0, downstream = 0)
{
      bed <- read.table(inputFilePath, sep="\t", header=header, stringsAsFactors = FALSE)
      if (dim(bed)[2] < 3)
      {
            stop("inputfile specified as ", inputFilePath, " is not a valid bed file!")
      }
      if (dim(bed)[2] >=6)
      {
           strand <- bed[,6]
      }
      else
      {
           strand <- rep("+", dim(bed)[1])
      }
      if (dim(bed)[2] >=4)
      {
           names <- bed[,4]
      }
      else
      {
           names <- paste(bed[,1], paste(bed[,2], bed[,3],sep="-"), sep=":")
      }
      chr <- bed[,1]
      Start <- as.numeric(bed[,2])
      End <- as.numeric(bed[,3])
      Start <- ifelse(strand=="-",
              Start - downstream,
              Start - upstream)
      End <- ifelse(strand=="-",
              End + upstream,
              End + downstream)
      starts <- unlist(apply(cbind(Start,1), 1, max))
      ends <- unlist(apply(cbind(End, seqlengths(BSgenomeName)[chr]), 1,min))
      seq <- BSgenome::getSeq(BSgenomeName, names = chr, start = starts,
            end = ends, strand = strand, width = NA, as.character = FALSE)
      names(seq) <- names
      seq
}
