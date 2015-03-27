uniqueREs <- function(REcutDetails, summary, offTargets, 
	scanUpstream = 100, scanDownstream = 100, BSgenomeName)
{
   	gRNAs <- as.character(summary$gRNAsPlusPAM)
	REs <- character(length(gRNAs))
	REwithName <- unique(cbind(as.character(REcutDetails$REpattern), 
		as.character(REcutDetails$REname)))
	onTargets <- subset(offTargets, as.character(offTargets$OffTargetSequence) %in% gRNAs)
	if (dim(onTargets)[1] >0)
	{
	   Start <- as.numeric(as.character(onTargets$chromStart)) - scanUpstream
	   End <- as.numeric(as.character(onTargets$chromEnd)) + scanDownstream
	   strand <- as.character(onTargets$strand)		
	   chr <- as.character(onTargets$chrom)
	   scanSequence = character(length(Start))
	   for (i in 1:length(Start))		
	   {
	      thisChr <-chr[i]
   	      thisEnd <- min(End[i], seqlengths(BSgenomeName)[thisChr][[1]])
	      thisStart <- max(1, Start[i])
 	      thisStrand <- as.character(strand[i])
	      scanSequence[i] <- getSeq(BSgenomeName, thisChr, start = thisStart, 
	         end = thisEnd, strand = thisStrand, 
	         as.character = TRUE)
	   }
	   for (i in 1:length(gRNAs))
	   {
		seq <- as.character(scanSequence[as.character(onTargets$OffTargetSequence) == as.character(gRNAs[i])])
		REnames <- unlist(strsplit(as.character(summary$REname[i]), " "))
		REpatterns <- unique(
	           REcutDetails$REpattern[as.character(REcutDetails$REname) %in% REnames])
		REnames <- REwithName[REwithName[,1] %in% REpatterns,2]
		REs[i] <- paste(unique(
		   REnames[isPatternUnique(seq, DNAStringSet(REpatterns)) == "Yes"]), collapse= " ")
	   }
	}
	else
	{
		REs = character(length(gRNAs))	
	}
	REs
}
