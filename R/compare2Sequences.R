compare2Sequences <- function(inputFile1Path, inputFile2Path, format = "fasta",
    findgRNAsWithREcutOnly = FALSE, REpatternFile, minREpatternSize = 6, 
    overlap.gRNA.positions = c(17, 18), findPairedgRNAOnly = FALSE, 
    min.gap = 0, max.gap = 20, gRNA.name.prefix = "gRNA", PAM.size = 3, 
    gRNA.size = 20, PAM = "NGG", PAM.pattern = "N[A|G]G$", max.mismatch = 3, 
    outputDir, weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 
    0.685, 0.583), overwrite = FALSE)
{
	append = ifelse(overwrite, FALSE, TRUE)
    gRNAs1 = offTargetAnalysis(inputFile1Path, format = format, 
        findgRNAs = TRUE,
        findPairedgRNAOnly = findPairedgRNAOnly, chromToSearch = "", 
        findgRNAsWithREcutOnly = findgRNAsWithREcutOnly, 
        REpatternFile = REpatternFile, minREpatternSize = minREpatternSize, 
        overlap.gRNA.positions =  overlap.gRNA.positions,
        min.gap = min.gap, max.gap = max.gap, 
        gRNA.name.prefix = gRNA.name.prefix, PAM.size = PAM.size,
        gRNA.size = gRNA.size, PAM = PAM, PAM.pattern = PAM.pattern,
        outputDir = file.path(outputDir, basename(inputFile1Path)), 
        weights = weights, overwrite = overwrite)
    gRNAs2 = offTargetAnalysis(inputFile2Path, format = format,                    
        findgRNAs = TRUE,
        findPairedgRNAOnly = findPairedgRNAOnly, chromToSearch = "",
        findgRNAsWithREcutOnly = findgRNAsWithREcutOnly, 
        REpatternFile = REpatternFile, minREpatternSize = minREpatternSize,
        overlap.gRNA.positions =  overlap.gRNA.positions, 
        min.gap = min.gap, max.gap = max.gap, 
        gRNA.name.prefix = gRNA.name.prefix, PAM.size = PAM.size,
        gRNA.size = gRNA.size, PAM = PAM, PAM.pattern = PAM.pattern, 
        outputDir = file.path(outputDir, basename(inputFile2Path)), 
            weights = weights, overwrite = overwrite)
	print("Scoring ...")
	subjects1 <- readDNAStringSet(inputFile1Path, format = "fasta",
        use.names = TRUE)
    subjects2 <- readDNAStringSet(inputFile2Path, format = "fasta",
        use.names = TRUE)
    outfile <- tempfile(tmpdir = getwd())
    max.mismatch <- max.mismatch + PAM.size -1
    revsubject <- reverseComplement(subjects2[[1]])
	chrom.len = nchar(as.character(subjects2))
	seqname <- names(subjects2)
    for (i in 1:length(gRNAs1))
    {
	    patternID <- gsub("'", "", names(gRNAs1)[i])
        if (length(patternID) < 1) {
           patternID <- paste("pattern", i, sep = "")
        }
        pattern <- DNAString(toupper(gRNAs1[[i]]))
        ### by default PAM is NGG or NAG
        plus_matches <- matchPattern(pattern, subjects2[[1]],
            max.mismatch = max.mismatch, min.mismatch = 0, 
            with.indels = FALSE, fixed = TRUE, algorithm = "auto")
        if (length(plus_matches) > 0) {
            names(plus_matches) <- rep.int(patternID, length(plus_matches))
            writeHits(pattern, seqname, plus_matches, strand = "+", 
                file = outfile, gRNA.size = gRNA.size,
                PAM = PAM.pattern, max.mismatch = max.mismatch - 2,
                chrom.len = chrom.len, append = append)
                append <- TRUE
		}
		if (reverseComplement(pattern) != pattern) {
            minus_matches <- matchPattern(pattern, revsubject, 
                max.mismatch = max.mismatch, min.mismatch = 0, 
                with.indels = FALSE, fixed = TRUE, algorithm = "auto")
            if (length(minus_matches) > 0) {
                names(minus_matches) <- rep.int(patternID,
                    length(minus_matches))
                writeHits(pattern, seqname, minus_matches, strand = "-", 
                    file = outfile, gRNA.size = gRNA.size, PAM = PAM.pattern,
                    max.mismatch = max.mismatch - 2, 
                    chrom.len = chrom.len, append = append)
                append <- TRUE
            }
        }
	} #### end of for loop for gRNAs1
   	seqname <- names(subjects1)
	revsubject <- reverseComplement(subjects1[[1]])
	chrom.len = nchar(as.character(subjects1))
    for (i in 1:length(gRNAs2))
    {
	    patternID <- gsub("'", "", names(gRNAs2)[i])
	    if (length(patternID) < 1) {
		    patternID <- paste("pattern", i, sep = "")
	    }
	    pattern <- DNAString(toupper(gRNAs2[[i]]))
        ### by default PAM is NGG or NAG
	    plus_matches <- matchPattern(pattern, subjects1[[1]],
			max.mismatch = max.mismatch, min.mismatch = 0, 
			with.indels = FALSE, fixed = TRUE, algorithm = "auto")
	    if (length(plus_matches) > 0) {
		    names(plus_matches) <- rep.int(patternID, length(plus_matches))
		    writeHits(pattern, seqname, plus_matches, strand = "+", 
			    file = outfile, gRNA.size = gRNA.size,
			    PAM = PAM.pattern, max.mismatch = max.mismatch - 2,
			    chrom.len = chrom.len, append = append)
		    append <- TRUE
	    }            
	    if (reverseComplement(pattern) != pattern) {
		    minus_matches <- matchPattern(pattern, revsubject, 
				max.mismatch = max.mismatch, min.mismatch = 0, 
				with.indels = FALSE, fixed = TRUE, algorithm = "auto")
		    if (length(minus_matches) > 0) {
			    names(minus_matches) <- rep.int(patternID,
				    length(minus_matches))
			    writeHits(pattern, seqname, minus_matches, strand = "-", 
				    file = outfile, gRNA.size = gRNA.size, PAM = PAM.pattern,
				    max.mismatch = max.mismatch - 2, 
				    chrom.len = chrom.len, append = append)
			    append <- TRUE
		    }
	    }
   } #### end of for loop for gRNAs2
    hits <- read.table(outfile, sep="\t", header = TRUE, 
        stringsAsFactors = FALSE)
	unlink(outfile)
	featureVectors <- buildFeatureVectorForScoring(hits = hits, 
		canonical.PAM = PAM, gRNA.size = gRNA.size)
	scores <- getOfftargetScore(featureVectors, weights = weights)
	targetInSeq1 <- scores$gRNAPlusPAM
	targetInSeq2 <- scores$gRNAPlusPAM
	scoreForSeq1 <- rep(100, dim(scores)[1])
	scoreForSeq2 <- rep(100, dim(scores)[1])
	
	scoreForSeq1[scores$chrom == names(subjects1)] <- 
		scores$score[scores$chrom == names(subjects1)]
	scoreForSeq2[scores$chrom == names(subjects2)] <- 
		scores$score[scores$chrom == names(subjects2)]
	targetInSeq1[scores$chrom == names(subjects1)] <- 
	    scores$OffTargetSequence[scores$chrom == names(subjects1)]
	targetInSeq2[scores$chrom == names(subjects2)] <- 
	    scores$OffTargetSequence[scores$chrom == names(subjects2)]
	seqNames <- c(names(subjects1), names(subjects2))
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
		targetSeqName = unlist(lapply(
		scores$chrom, function(i) {seqNames[seqNames !=i]}))
	)
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
			targetSeqName = names(subjects1)
			)
		seqs.new <- rbind(seqs.new, seqs1.only)
	}
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
			targetSeqName = names(subjects2)
		)
		seqs.new <- rbind(seqs.new, seqs2.only)
	}
	seqs = cbind(seqs.new, 
        scoreDiff = round(as.numeric(seqs.new[,7]) - as.numeric(seqs.new[,8]),4))
    setwd(outputDir)
    write.table(seqs[order(as.numeric(seqs[,12]), decreasing = TRUE), ], 
        file = "scoresFor2InputSequences.xls",
        sep = "\t", row.names = FALSE, col.names=TRUE)
	print("Done!")
    seqs
}
