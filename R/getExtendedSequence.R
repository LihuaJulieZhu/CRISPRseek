getExtendedSequence <- function(targets, BSgenomeName, baseBeforegRNA = 13, baseAfterPAM = 24, forMethod = "Lindel")
{
   if (forMethod == "Lindel")
   {
        baseBeforegRNA = 13
        baseAfterPAM = 24
   }
   if (dim(targets)[1] > 0)
   {
        chr <- as.character(targets$chrom)
        strand <- as.character(targets$strand)
        Start <- ifelse(strand=="-",
              as.numeric(as.character( targets$chromStart)) - baseAfterPAM,
              as.numeric(as.character( targets$chromStart)) - baseBeforegRNA)
        End <- ifelse(strand=="-",
              as.numeric(as.character( targets$chromEnd)) + as.numeric(baseBeforegRNA),
              as.numeric(as.character( targets$chromEnd)) + as.numeric(baseAfterPAM))
    }
    starts <- unlist(apply(cbind(Start,1), 1, max))
    ends <- unlist(apply(cbind(End, seqlengths(BSgenomeName)[chr]), 1,min))
        
    extendedSequence <- getSeq(BSgenomeName, names = chr, start = starts,
           end = ends, strand = strand, width = NA, as.character = TRUE)
}
