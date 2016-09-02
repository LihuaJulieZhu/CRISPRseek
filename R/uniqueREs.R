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
	}
	REs
}
