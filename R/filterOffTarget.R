filterOffTarget <-
    function(scores, min.score = 0.5, topN = 100, topN.OfftargetTotalScore = 10,
    annotateExon = TRUE, txdb, orgAnn, outputDir, oneFilePergRNA = FALSE,
    fetchSequence = TRUE, upstream = 200, downstream = 200, BSgenomeName,
	baseBeforegRNA = 4, baseAfterPAM = 3,
	featureWeightMatrixFile = system.file("extdata", "DoenchNBT2014.csv", 
		package = "CRISPRseek"))
{
    if (featureWeightMatrixFile != ""  & file.exists(featureWeightMatrixFile))
	{
		featureWeightMatrix <- read.csv(featureWeightMatrixFile, header=TRUE)
	}
	if (fetchSequence && (missing(BSgenomeName) || 
        class(BSgenomeName) != "BSgenome")) 
    {
        stop("To fetch sequences, BSgenomeName is required as BSgenome object!")
    }
    if (annotateExon && ( missing(txdb) || (class(txdb) != "TxDb" && class(txdb) != "TranscriptDb")))
    {
        stop("To indicate whether an offtarget is inside an exon, txdb is
            required as TxDb object!")
    }
    if (annotateExon && missing(orgAnn))
    {
	    warning("orgAnn was not included. See the updated manual for information about how to use the orgAnn parameter to generate gene identifiers in the offTarget output file")
    }
    scores <- scores[scores$score >= min.score,]
    if (length(grep("IsMismatch.pos", colnames(scores))) > 0)
        scores <- scores[,-c(grep("IsMismatch.pos", colnames(scores)))]
	if (substr(outputDir, nchar(outputDir), nchar(outputDir)) != .Platform$file.sep)
    {
        outputDir <- paste(outputDir, "", sep = .Platform$file.sep)
    }	
    if ( ! file.exists(outputDir))
    {
        dir.create(outputDir)
    }
    OfftargetFile <-paste(outputDir, "OfftargetAnalysis.xls", sep = "")
    OfftargetSummary <-paste(outputDir, "Summary.xls", sep = "")
    gRNAsPlusPAM<- unique(cbind(gRNAPlusPAM = scores$gRNAPlusPAM, 
	forViewInUCSC =scores$forViewInUCSC))
    names <- gRNAsPlusPAM[,1]
    top5OfftargetTotalScore <- numeric(length(names))
    topNOfftargetTotalScore <- top5OfftargetTotalScore
    temp <- cbind(names, forViewInUCSC = gRNAsPlusPAM[,2], top5OfftargetTotalScore, 
        topNOfftargetTotalScore)
    mismatch.distance2PAM <- matrix(ncol = 11, nrow = length(names))
    append <- FALSE
    for (i in 1:length(names))
    {
        this.score <- scores[scores$gRNAPlusPAM == names[i],]
        this.score <- this.score[order(this.score$score, decreasing = TRUE),]
        maxN <- min(topN+1, dim(this.score)[1])
        this.score <- this.score[1:maxN,]
        maxN.totalScore <- min(maxN, (topN.OfftargetTotalScore + 1))
        if (maxN < 6)
            temp[i,3] <- sum(this.score$score[2:maxN])
        else
            temp[i,3] <- sum(this.score$score[2:6]) 
	    if (maxN < maxN.totalScore)
	        temp[i,4] <- sum(this.score$score[2:maxN])
	    else
            temp[i,4] <- sum(this.score$score[2:maxN.totalScore])
        temp[i,1] <- unique(this.score$name)
	    forSummary <- this.score[1:11,]
	    forSummary <- forSummary[order(forSummary$n.mismatch),]
        mismatch.distance2PAM[i,] <- 
            as.character(forSummary$mismatch.distance2PAM)
	    forSummary <- forSummary[2:11,]
	    forSummary <- forSummary[order(forSummary$score, decreasing=TRUE),]
	    mismatch.distance2PAM[i,2:11] <-
            as.character(forSummary$mismatch.distance2PAM)
		this.score <- cbind(name = this.score$name, 
                gRNAPlusPAM = this.score$gRNAPlusPAM,
                OffTargetSequence = this.score$OffTargetSequence,
                score = this.score$score, n.mismatch = this.score$n.mismatch, 
                mismatch.distance2PAM = 
                as.character(this.score$mismatch.distance2PAM), 
                alignment = this.score$alignment,
                NGG = as.character(this.score$NGG),
                forViewInUCSC = this.score$forViewInUCSC, 
                strand = this.score$strand, 
                chrom = this.score$chrom, chromStart = this.score$chromStart,
                chromEnd = this.score$chromEnd)
        if (oneFilePergRNA & dim(this.score)[1] > 0)
            write.table(this.score[!is.na(this.score[,grep("score", 
                colnames(this.score))]),], 
                file = paste( outputDir, "OfftargetAnalysis-", 
                as.character(temp[i,1]), ".xls", sep = ""), sep = "\t",
                row.names = FALSE)
        if (i == 1 && dim(this.score)[1] > 0)
        {
            write.table(this.score, file = OfftargetFile, sep = "\t", 
                row.names = FALSE,  append = append)
            append <- TRUE
        }
        else if (dim(this.score)[1] > 0)
        {
            #this.score <- this.score[!is.na(this.score[,grep("score", 
            ##colnames(this.score))]),]
            write.table(this.score, file = OfftargetFile, sep = "\t",
                row.names = FALSE, col.names = FALSE, append = append)
            append <- TRUE
        }
    }
    temp <- cbind(temp, mismatch.distance2PAM)
    colnames(temp)[5] <- "top1Hit(onTarget)MMdistance2PAM"
    colnames(temp)[4] <- paste("top", topN.OfftargetTotalScore, 
       "OfftargetTotalScore", sep = "")
    colnames(temp)[6:15] <- paste("topOfftarget", 1:10, "MMdistance2PAM",
        sep = "")
    Offtargets <- read.table(OfftargetFile, sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE)
	if (annotateExon)
	{
		Offtargets <- annotateOffTargets(Offtargets, txdb, orgAnn)
	}
	Start <- as.numeric(as.character(Offtargets$chromStart))
	End <- as.numeric(as.character(Offtargets$chromEnd))
	strand <- as.character(Offtargets$strand)		
	Start[strand == "+"] = Start[strand == "+"] -  baseBeforegRNA
	Start[strand == "-"] = Start[strand == "-"] -  baseAfterPAM
	End[strand == "+"] = End[strand == "+"] +  baseAfterPAM
        End[strand == "-"] = End[strand == "-"] +  baseBeforegRNA
       
	chr <- as.character(Offtargets$chrom)
	for (i in 1:length(Start))		
	{
		thisChr <-chr[i]
		thisEnd <- min(End[i], seqlengths(BSgenomeName)[thisChr][[1]])
		thisStart <- max(1, Start[i])
		thisStrand <- as.character(strand[i])
		if (i == 1)
		{
			extendedSequence <- getSeq(BSgenomeName, thisChr, start = thisStart, 
			   end = thisEnd, strand = thisStrand, width = NA, 
			   as.character = TRUE)
		}
		else
		{
			extendedSequence <- c(extendedSequence, 
				 getSeq(BSgenomeName, thisChr, start = thisStart, 
	 			 end = thisEnd, strand = thisStrand, width = NA, 
				 as.character = TRUE))
		}
	}
	Offtargets <- cbind(Offtargets, extendedSequence = extendedSequence)
	gRNAefficiency <- calculategRNAEfficiency(extendedSequence, 
		baseBeforegRNA = baseBeforegRNA,
		featureWeightMatrix = featureWeightMatrix)
	Offtargets <- cbind(Offtargets, gRNAefficacy = gRNAefficiency)
	if (fetchSequence)
	{
		Start <- as.numeric(as.character(Offtargets$chromStart)) - as.numeric(upstream)
        End <- as.numeric(as.character(Offtargets$chromEnd)) + as.numeric(downstream)	
        strand <- as.character(Offtargets$strand)		
        chr <- as.character(Offtargets$chrom)
        for (i in 1:length(Start))		
	   {
		   thisChr <-chr[i]
		   thisEnd <- min(End[i], seqlengths(BSgenomeName)[thisChr][[1]])
		   thisStart <- max(1, Start[i])
		   thisStrand <- as.character(strand[i])
		   if (i == 1)
		   {
			   seq <- getSeq(BSgenomeName, thisChr, start = thisStart, 
				  end = thisEnd, strand = thisStrand, width = NA, 
				  as.character = TRUE)
		   }
		   else
		   {
			  seq <- c(seq, getSeq(BSgenomeName, thisChr, start = thisStart, 
				  end = thisEnd, strand = thisStrand, width = NA, 
				  as.character = TRUE))
		   }
	   }
	   Offtargets <- cbind(Offtargets, flankSequence = seq)
	}
	
    write.table(temp, file = OfftargetSummary, sep = "\t", row.names = FALSE)
    write.table(unique(Offtargets[order(as.character(Offtargets$name), 
        -as.numeric(as.character(Offtargets$score)), 
        as.character(Offtargets$OffTargetSequence)),]), 
        file = OfftargetFile, sep = "\t", row.names = FALSE)
    list(offtargets = unique(Offtargets), summary = temp)
}
