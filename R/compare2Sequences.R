compare2Sequences <- function(inputFile1Path, inputFile2Path, format = "fasta",
    findSpacersWithREcutOnly = FALSE, REpatternFile, minREpatternSize = 6, 
    overlap.spacer.positions = c(17, 18), findPairedSpacerOnly = FALSE, 
    min.gap = 0, max.gap = 20, spacer.name.prefix = "gRNA", PAM.size = 3, 
    spacer.size = 20, PAM = "NGG", PAM.pattern = "N[A|G]G$", outputDir, 
    weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 
    0.685, 0.583), overwrite = FALSE)
{
    weights = rev(weights)
    spacers1 = offTargetAnalysis(inputFile1Path, format = format, 
        findSpacers = TRUE,
        findPairedSpacerOnly = findPairedSpacerOnly, chromToSearch = "", 
        findSpacersWithREcutOnly = findSpacersWithREcutOnly, 
        REpatternFile = REpatternFile, minREpatternSize = minREpatternSize, 
        overlap.spacer.positions =  overlap.spacer.positions,
        min.gap = min.gap, max.gap = max.gap, 
        spacer.name.prefix = spacer.name.prefix, PAM.size = PAM.size,
        spacer.size = spacer.size, PAM = PAM, PAM.pattern = PAM.pattern,
        outputDir = outputDir, weights = weights, overwrite = overwrite)
    spacers2 = offTargetAnalysis(inputFile2Path, format = format,                    
        findSpacers = TRUE,
        findPairedSpacerOnly = findPairedSpacerOnly, chromToSearch = "",
        findSpacersWithREcutOnly = findSpacersWithREcutOnly, 
        REpatternFile = REpatternFile, minREpatternSize = minREpatternSize,
        overlap.spacer.positions =  overlap.spacer.positions, 
        min.gap = min.gap, max.gap = max.gap, 
        spacer.name.prefix = spacer.name.prefix, PAM.size = PAM.size,
        spacer.size = spacer.size, PAM = PAM, PAM.pattern = PAM.pattern, 
        outputDir = outputDir, weights = weights, overwrite = overwrite)
    seqs1 = cbind(names(spacers1), as.character(spacers1))
    seqs2 = cbind(names(spacers2), as.character(spacers2))

    colnames(seqs1) = c("name", "spacerPlusPAM")
    colnames(seqs2) = c("name", "spacerPlusPAM")
    seqs = merge(seqs1, seqs2, by = "spacerPlusPAM", all = TRUE)
    name = as.character(seqs[,2])
    name[is.na(name)] = as.character(seqs[is.na(seqs[,2]),3])
    seqs = cbind(name, seqs)
    colnames(seqs)[2] = "spacerPlusPAM"
    seqs[,3:4] = cbind(targetInSeq1 = character(dim(seqs)[1]), 
        targetInSeq2 = character(dim(seqs)[1]))
    colnames(seqs)[3:4] = c("targetInSeq1", "targetInSeq2")
    subjects1 <- readDNAStringSet(inputFile1Path, format = "fasta",
        use.names = TRUE)
    subjects2 <- readDNAStringSet(inputFile2Path, format = "fasta",
        use.names = TRUE)
    seqs = cbind(seqs, scoreForSeq1 = numeric(dim(seqs)[1]), 
        scoreForSeq2 = numeric(dim(seqs)[1]))
    seqs = cbind(seqs, mismatch.distance2PAM = numeric(dim(seqs)[1]), 
         n.mismatch = numeric(dim(seqs)[1]))
    #IsMismatch.posX = matrix(0, ncol = spacer.size, nrow = (dim(seqs)[1]) 
    targetSeqName = character(dim(seqs)[1]) 
    for (i in 1:length(spacers1))
    {
        if (length(grep(paste(spacer.name.prefix,"f", sep = ""),
            names(spacers1)[i])))
	    aln <- pairwiseAlignment(spacers1[i], subjects2[1], 
                type = "global-local")
	else
	    aln <- pairwiseAlignment(spacers1[i], 
                reverseComplement(subjects2[1]), type = "global-local")
	if (length(mismatchTable(aln)$PatternStart))
       	    seqs$mismatch.distance2PAM[i] <- width(spacers1[i]) -
                mismatchTable(aln)$PatternStart - PAM.size  + 1
	else
	    seqs$mismatch.distance2PAM[i] <- 0
	seqs$targetInSeq2[i] <- as.character(subject(aln))
	seqs$targetInSeq1[i] <- as.character(pattern(aln))
	seqs$name[i] <- names(spacers1)[i]
	seqs$spacerPlusPAM[i] <- as.character(spacers1)[i]
	seqs$scoreForSeq1[i] <- 100
	seqs$n.mismatch[i] = nedit(aln)
	if (seqs$mismatch.distance2PAM[i] >0)
	    seqs$scoreForSeq2[i] <- 
                100 * (1 - weights[seqs$mismatch.distance2PAM[i]])
	else
	    seqs$scoreForSeq2[i] <- 100
    }

    ioffset <- length(spacers1)
    targetSeqName[1:ioffset] <- names(subjects1)
    for (i in 1:length(spacers2))
    {
	j <- i + ioffset
	if (length(grep(paste(spacer.name.prefix,"f", sep = ""), 
            names(spacers2)[i])))
            aln <- pairwiseAlignment(spacers2[i], subjects1[1], 
                type = "global-local")
	else
	    aln <- pairwiseAlignment(spacers2[i], 
                reverseComplement(subjects1[1]),
	        type = "global-local")
	if (length(mismatchTable(aln)$PatternStart))
       	    seqs$mismatch.distance2PAM[j] <- width(spacers2[i]) - 
                mismatchTable(aln)$PatternStart - PAM.size  + 1
	else
            seqs$mismatch.distance2PAM[j] <- 0 
	seqs$targetInSeq2[j] <- as.character(subject(aln))
	seqs$targetInSeq1[j] <- as.character(pattern(aln))
	seqs$name[j] <- names(spacers2)[i]
	seqs$spacerPlusPAM[j] <- as.character(spacers2)[i]
 	seqs$scoreForSeq2[j] <- 100
	seqs$n.mismatch[j] = nedit(aln)
	if (seqs$mismatch.distance2PAM[j] >0)
	    seqs$scoreForSeq1[j] <- 
                100 * (1 - weights[seqs$mismatch.distance2PAM[j]])
	else
	     seqs$scoreForSeq1[j] <- 100
    }
    j <- ioffset + 1
    k <- dim(seqs)[1]
    targetSeqName[j:k] <- names(subjects2)
    seqs = cbind(seqs, scoreDiff = seqs$scoreForSeq1 - seqs$scoreForSeq2, targetSeqName)
    setwd(outputDir)
    write.table(seqs[order(seqs$scoreDiff, decreasing = TRUE), ], 
        file = "scoresFor2InputSequences.xls",
        sep = "\t", row.names = FALSE)
    seqs
}
