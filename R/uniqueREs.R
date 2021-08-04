#' Output restriction enzymes that recognize only the gRNA cleavage sites
#' 
#' For each identified gRNA, output restriction enzymes that recognize only the
#' gRNA cleavage sites.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param REcutDetails REcutDetails stored in the REcutDetails.xls
#' @param summary summary stored in the summary.xls
#' @param offTargets offTargets stored in the offTargets.xls
#' @param scanUpstream upstream offset from the gRNA start, default 100
#' @param scanDownstream downstream offset from the gRNA end, default 100
#' @param BSgenomeName BSgenome object. Please refer to available.genomes in
#' BSgenome package. For example, 
#' \itemize{
#' \item{BSgenome.Hsapiens.UCSC.hg19} - {for hg19}
#' \item{BSgenome.Mmusculus.UCSC.mm10} - {for mm10}
#' \item{BSgenome.Celegans.UCSC.ce6} - {for ce6}
#' \item{BSgenome.Rnorvegicus.UCSC.rn5} - {for rn5}
#' \item{BSgenome.Drerio.UCSC.danRer7} - {for Zv9}
#' \item{BSgenome.Dmelanogaster.UCSC.dm3} - {for dm3}
#' }
#' @return returns the RE sites that recognize only the gRNA cleavage sites for
#' each gRNA.
#' @note %% ~~further notes~~
#' @author Lihua Julie Zhu
#' @seealso
#' @references %% ~put references to the literature/web site here ~
#' @keywords misc
#' @examples
#' 
#'     library("BSgenome.Hsapiens.UCSC.hg19")
#'     load(system.file("extdata", "ForTestinguniqueREs.RData",
#'             package = "CRISPRseek"))
#'     uniqueREs(results$REcutDetails, results$summary, results$offtarget,
#' 	scanUpstream = 50,
#'         scanDownstream = 50, BSgenomeName = Hsapiens)
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom Biostrings DNAStringSet
#' @importFrom BiocGenerics cbind unlist
#' @export
uniqueREs <- function(REcutDetails, summary, offTargets, 
	scanUpstream = 100, scanDownstream = 100, BSgenomeName)
{
	REwithName <- unique(cbind(as.character(REcutDetails$REpattern), 
		as.character(REcutDetails$REname)))
        REs = character(dim(summary)[1])
        summary$id <- summary$names
        summary$id <- paste(summary$names, summary$forViewInUCSC, sep="-")
        offTargets$id <- offTargets$name
        offTargets$id <- paste(offTargets$name, offTargets$forViewInUCSC, sep="-")
        summary <- cbind(summary, 
            offTargets[match(summary$id, offTargets$id),
                c("chrom", "chromStart", "chromEnd", "strand")])
        REs = character(dim(summary)[1])
	if (dim(summary)[1] >0)
	{
	   Start <- as.numeric(as.character(summary$chromStart)) - scanUpstream
	   End <- as.numeric(as.character(summary$chromEnd)) + scanDownstream
	   strand <- as.character(summary$strand)		
	   chr <- as.character(summary$chrom)
	   for (i in 1:length(Start))		
	   {
	       thisChr <-chr[i]
               if (!is.na(thisChr) && thisChr != "")
               {
   	          thisEnd <- min(End[i], seqlengths(BSgenomeName)[thisChr][[1]])
	          thisStart <- max(1, Start[i])
 	          thisStrand <- as.character(strand[i])
	          scanSequence <- BSgenome::getSeq(BSgenomeName, 
                     thisChr, start = thisStart, 
	             end = thisEnd, strand = thisStrand, 
	             as.character = TRUE)
                  REnames <- unlist(strsplit(as.character(summary$REname[i]), " "))
	          REpatterns <- unique(
	              REcutDetails$REpattern[as.character(
                       REcutDetails$REname) %in% REnames])
	          REnames <- REwithName[REwithName[,1] %in% REpatterns,2]
	          REs[i] <- paste(unique(
		     REnames[isPatternUnique(scanSequence, 
                     DNAStringSet(REpatterns)) == "Yes"]), collapse= " ")
               }
               else
               {
                 REs[i] = ""
               }
	   }
	}
	REs
}
