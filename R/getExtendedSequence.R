#' @importFrom BiocGenerics unlist cbind
#' @importFrom BSgenome getSeq
#' @importFrom IRanges width
#' @importFrom methods as
#' @importFrom GenomeInfoDb seqlengths
getExtendedSequence <- function(targets, useBSgenome = TRUE, BSgenomeName, 
     genomeSeq, baseBeforegRNA = 13, baseAfterPAM = 24, forMethod = "Lindel")
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
    if (useBSgenome)
    {
         ends <- unlist(apply(cbind(End, seqlengths(BSgenomeName)[chr]), 1,min))
         extendedSequence <- getSeq(BSgenomeName, names = chr, start = starts,
             end = ends, strand = strand, width = NA, as.character = TRUE)
     }
    else
    {
         ends <- unlist(apply(cbind(End, width(genomeSeq)[names(genomeSeq) %in% chr]),
                1, min))
         extended.info <- data.frame(chrom = chr, start = starts, end = ends, strand = strand)
         extendedSequence <- as.character(getSeq(genomeSeq, as(extended.info, "GRanges")))
    }
    extendedSequence
}
