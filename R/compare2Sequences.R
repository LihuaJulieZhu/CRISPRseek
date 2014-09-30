compare2Sequences <- function(inputFile1Path, inputFile2Path, format = "fasta",
    findgRNAsWithREcutOnly = FALSE, searchDirection=c("both","1to2", "2to1"),
    REpatternFile=system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek"),
    minREpatternSize = 6, 
    overlap.gRNA.positions = c(17, 18), findPairedgRNAOnly = FALSE, 
    min.gap = 0, max.gap = 20, gRNA.name.prefix = "gRNA", PAM.size = 3, 
    gRNA.size = 20, PAM = "NGG", PAM.pattern = "N[A|G]G$", max.mismatch = 3, 
    outputDir, weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 
    0.685, 0.583), overwrite = FALSE)
{
	append = ifelse(overwrite, FALSE, TRUE)
	if (class(inputFile1Path) != "DNAStringSet")
		outputDir1 <- file.path(outputDir, basename(inputFile1Path))
	else
		outputDir1 <- file.path(outputDir, "File1")
	if (class(inputFile2Path) != "DNAStringSet")
		outputDir2 <- file.path(outputDir, basename(inputFile2Path))
	else
		outputDir2 <- file.path(outputDir, "File2")
	if(searchDirection == "both" || searchDirection == "1to2")
	{
    gRNAs1 = offTargetAnalysis(inputFile1Path, format = format, 
        findgRNAs = TRUE,
        findPairedgRNAOnly = findPairedgRNAOnly, chromToSearch = "", 
        findgRNAsWithREcutOnly = findgRNAsWithREcutOnly, 
        REpatternFile = REpatternFile, minREpatternSize = minREpatternSize, 
        overlap.gRNA.positions =  overlap.gRNA.positions,
        min.gap = min.gap, max.gap = max.gap, 
        gRNA.name.prefix = gRNA.name.prefix, PAM.size = PAM.size,
        gRNA.size = gRNA.size, PAM = PAM, PAM.pattern = PAM.pattern,
        outputDir = outputDir1, 
        weights = weights, overwrite = overwrite)
	}
	if(searchDirection == "both" || searchDirection == "2to1")
	{
    gRNAs2 = offTargetAnalysis(inputFile2Path, format = format,                    
        findgRNAs = TRUE,
        findPairedgRNAOnly = findPairedgRNAOnly, chromToSearch = "",
        findgRNAsWithREcutOnly = findgRNAsWithREcutOnly, 
        REpatternFile = REpatternFile, minREpatternSize = minREpatternSize,
        overlap.gRNA.positions =  overlap.gRNA.positions, 
        min.gap = min.gap, max.gap = max.gap, 
        gRNA.name.prefix = gRNA.name.prefix, PAM.size = PAM.size,
        gRNA.size = gRNA.size, PAM = PAM, PAM.pattern = PAM.pattern, 
        outputDir = outputDir2, 
            weights = weights, overwrite = overwrite)
	}
    print("Scoring ...")
    if (class(inputFile1Path) != "DNAStringSet")
    {
    	subjects1 <- readDNAStringSet(inputFile1Path, format = format,
    	    use.names = TRUE)
    }
    else
    {
	subjects1 <- inputFile1Path
    }
   if (class(inputFile2Path) != "DNAStringSet")
   {
    	subjects2 <- readDNAStringSet(inputFile2Path, format = format,
       	    use.names = TRUE)
   }
   else
   {
	subjects2 <- inputFile2Path
   }
    outfile <- tempfile(tmpdir = getwd())
    max.mismatch <- max.mismatch + PAM.size -1
	seqname <- names(subjects2)
    seqname <- gsub("'", "", seqname)
    seqname <- gsub(" ", "", seqname)
    seqname <- gsub("\t", ":", seqname)
	names(subjects2) <- seqname
#revsubject <- reverseComplement(subjects2[[1]])
	revsubject <- reverseComplement(subjects2)
	chrom.len = nchar(as.character(subjects2))
    if(searchDirection == "both" || searchDirection == "1to2")
	{
    for (i in 1:length(gRNAs1))
    {
	    patternID <- gsub("'", "", names(gRNAs1)[i])
	    patternID <- gsub(" ", "", patternID)
	    patternID <- gsub("\t", ":", patternID)
        if (length(patternID) < 1) {
           patternID <- paste("pattern", i, sep = "")
        }
        pattern <- DNAString(toupper(gRNAs1[[i]]))
        ### by default PAM is NGG or NAG
		for (j in 1:length(subjects2))
		{
			plus_matches <- matchPattern(pattern, subjects2[[j]],
#plus_matches <- vmatchPattern(pattern, subjects2,
				max.mismatch = max.mismatch, min.mismatch = 0, 
                with.indels = FALSE, fixed = TRUE, algorithm = "auto")
			if (length(plus_matches) > 0) {
				names(plus_matches) <- rep.int(patternID, length(plus_matches))
				writeHits(pattern, seqname[j], plus_matches, strand = "+", 
					file = outfile, gRNA.size = gRNA.size,
					PAM = PAM.pattern, max.mismatch = max.mismatch - 2,
					chrom.len = as.numeric(chrom.len[j]), append = append)
				append <- TRUE
				plus_matches <- ""
			}
		}
		if (reverseComplement(pattern) != pattern) {
			for (j in 1:length(subjects2))
			{
				minus_matches <- matchPattern(pattern, revsubject[[j]], 
					max.mismatch = max.mismatch, min.mismatch = 0, 
					with.indels = FALSE, fixed = TRUE, algorithm = "auto")
				if (length(minus_matches) > 0) {
					names(minus_matches) <- rep.int(patternID,
						length(minus_matches))
					writeHits(pattern, seqname[j], minus_matches, strand = "-", 
						file = outfile, gRNA.size = gRNA.size, PAM = PAM.pattern,
						max.mismatch = max.mismatch - 2, 
						chrom.len = as.numeric(chrom.len[j]), append = append)
					append <- TRUE
					minus_mathes <- ""
				}
            }
        }
	} #### end of for loop for gRNAs1
	} # end of if searchDirection == "both" or searchDirection == "1to2"
	cat("finish off-target search in sequence 2\n") 
   	seqname <- names(subjects1)
	seqname <- gsub("'", "", seqname)
	seqname <- gsub(" ", "", seqname)
	seqname <- gsub("\t", ":", seqname)
	revsubject <- reverseComplement(subjects1)
	#revsubject <- reverseComplement(subjects1[[1]])
	names(subjects1) <- seqname

	chrom.len = nchar(as.character(subjects1))
	if(searchDirection == "both" || searchDirection == "2to1")
	{
    for (i in 1:length(gRNAs2))
    {
	    patternID <- gsub("'", "", names(gRNAs2)[i])
            patternID <- gsub(" ", "", patternID)
            patternID <- gsub("\t", ":", patternID)
	    if (length(patternID) < 1) {
		    patternID <- paste("pattern", i, sep = "")
	    }
	    pattern <- DNAString(toupper(gRNAs2[[i]]))
        ### by default PAM is NGG or NAG
		for (j in 1:length(subjects1))
		{
			plus_matches <- matchPattern(pattern, subjects1[[j]],
				max.mismatch = max.mismatch, min.mismatch = 0, 
				with.indels = FALSE, fixed = TRUE, algorithm = "auto")
			if (length(plus_matches) > 0) {
				names(plus_matches) <- rep.int(patternID, length(plus_matches))
				writeHits(pattern, seqname[j], plus_matches, strand = "+", 
					file = outfile, gRNA.size = gRNA.size,
					PAM = PAM.pattern, max.mismatch = max.mismatch - 2,
					chrom.len = as.numeric(chrom.len[j]), append = append)
				append <- TRUE
				plus_matches <- ""
			}
	    }            
	    if (reverseComplement(pattern) != pattern) {
			for (j in 1:length(subjects1))
			{				
				minus_matches <- matchPattern(pattern, revsubject[[j]], 
					max.mismatch = max.mismatch, min.mismatch = 0, 
					with.indels = FALSE, fixed = TRUE, algorithm = "auto")
				if (length(minus_matches) > 0) {
					names(minus_matches) <- rep.int(patternID,
				    length(minus_matches))
					writeHits(pattern, seqname[j], minus_matches, strand = "-", 
						file = outfile, gRNA.size = gRNA.size, PAM = PAM.pattern,
						max.mismatch = max.mismatch - 2, 
						chrom.len = as.numeric(chrom.len[j]), append = append)
					append <- TRUE
					minus_matches <- ""
				}
		    }
	    }
   } #### end of for loop for gRNAs2
	} # if searchDirection == "both" or searchDirection == "1to2"
	cat("finish off-target search in sequence 1\n")
	if (file.exists(outfile))
	{
		hits <- read.table(outfile, sep="\t", header = TRUE, 
        stringsAsFactors = FALSE)
#cat(dim(hits))
	}
	else
	{
		stop("No offtarget found! You can alter your 
			 search criteria such as increasing max.mismatch!")
	}
	unlink(outfile)
	featureVectors <- buildFeatureVectorForScoring(hits = hits, 
		canonical.PAM = PAM, gRNA.size = gRNA.size)
	cat("finish feature vector building\n")
	scores <- getOfftargetScore(featureVectors, weights = weights)
	cat("finish score calculation\n")
	targetInSeq1 <- scores$gRNAPlusPAM
	targetInSeq2 <- scores$gRNAPlusPAM
	scoreForSeq1 <- rep(100, dim(scores)[1])
	scoreForSeq2 <- rep(100, dim(scores)[1])
	
	if(length(subjects1) == 1 && length(subjects2) == 1)
	{
		scoreForSeq1[scores$chrom == names(subjects1)] <- 
			scores$score[scores$chrom == names(subjects1)]
		scoreForSeq2[scores$chrom == names(subjects2)] <- 
			scores$score[scores$chrom == names(subjects2)]
		targetInSeq1[scores$chrom == names(subjects1)] <- 
			scores$OffTargetSequence[scores$chrom == names(subjects1)]
		targetInSeq2[scores$chrom == names(subjects2)] <- 
			scores$OffTargetSequence[scores$chrom == names(subjects2)]
		seqNames <- c(names(subjects1), names(subjects2))
		targetSeqName = unlist(lapply(
			scores$chrom, function(i) {seqNames[seqNames !=i]}))
	}
	if (searchDirection == "1to2")
	{
		scoreForSeq2 = scores$score
		targetInSeq2 = scores$OffTargetSequence
		targetSeqName = names(subjects1)
	}
	if (searchDirection == "2to1")
	{
		scoreForSeq1 = scores$score
		targetInSeq1 = scores$OffTargetSequence
		targetSeqName = names(subjects2)
	}
	if(searchDirection == "both" && length(subjects1) > 1 && length(subjects2) > 1)
	{
		cat("Needs to be implemented yet")
	}
	
	seqs.new <- cbind(name = scores$name,
		gRNAPlusPAM = scores$gRNAPlusPAM,
		targetInSeq1 = targetInSeq1,
		targetInSeq2 = targetInSeq2,
	    guideAlignment2OffTarget = scores$alignment,
	    offTargetStrand = scores$strand,
		scoreForSeq1 = scoreForSeq1,
		scoreForSeq2 = scoreForSeq2,
		mismatch.distance2PAM = as.character(scores$mismatche.distance2PAM),
		n.mismatch = scores$n.mismatch,
		offTarget = scores$forViewInUCSC,
		targetSeqName = targetSeqName
	)
	if(searchDirection == "both" || searchDirection == "1to2")
	{
	if (length(setdiff(names(gRNAs1),seqs.new[,1])) >0)
	{
		seqs1.only <- cbind(name = names(gRNAs1)[!names(gRNAs1) %in% seqs.new[,1]],
			gRNAPlusPAM = as.character(gRNAs1)[!as.character(gRNAs1) %in% seqs.new[,2]],
			targetInSeq1 = as.character(gRNAs1)[!as.character(gRNAs1) %in% seqs.new[,2]],
			targetInSeq2 = "NA",
		    guideAlignment2OffTarget = "NA",
		    offTargetStrand = "NA",
			scoreForSeq1 = 100,
			scoreForSeq2 = 0,
			mismatch.distance2PAM = "NA",
			n.mismatch = "NA",
			offTarget = "NA",
			targetSeqName = names(subjects1)
			)
		seqs.new <- rbind(seqs.new, seqs1.only)
	}
	}
	if(searchDirection == "both" || searchDirection == "2to1")
	{
	if (length(setdiff(names(gRNAs2),seqs.new[,1])) >0)
	{
		seqs2.only <- cbind(name = names(gRNAs2)[!names(gRNAs2) %in% seqs.new[,1]],
    		gRNAPlusPAM = as.character(gRNAs2)[!as.character(gRNAs2) %in% seqs.new[,2]],
		    targetInSeq1 = "NA",
			targetInSeq2 = as.character(gRNAs2)[!as.character(gRNAs2) %in% seqs.new[,2]],
		    guideAlignment2OffTarget = "NA",
            offTargetStrand = "NA",
            scoreForSeq1 = 0,
			scoreForSeq2 = 100,
			mismatch.distance2PAM = "NA",
			n.mismatch = "NA",
			offTarget = "NA",
			targetSeqName = names(subjects2)
		)
		seqs.new <- rbind(seqs.new, seqs2.only)
	}
	}
	seqs = unique(cbind(seqs.new, 
        scoreDiff = round(as.numeric(seqs.new[,7]) - as.numeric(seqs.new[,8]),4)))
	
    setwd(outputDir)
	if (dim(seqs)[1] ==1)
	{
		write.table(seqs, file = "scoresFor2InputSequences.xls",
				sep = "\t", row.names = FALSE, col.names=TRUE)
	}
	else
	{
		write.table(seqs[order(as.numeric(seqs[,13]), decreasing = TRUE), ], 
			file = "scoresFor2InputSequences.xls",
			sep = "\t", row.names = FALSE, col.names=TRUE)
	}
	print("Done!")
#scores
 seqs
}
