filterOffTarget <-
    function(scores, min.score = 0.5, topN = 100, topN.OfftargetTotalScore = 10,
    annotateExon = TRUE, txdb, orgAnn, outputDir, oneFilePergRNA = FALSE,
    fetchSequence = TRUE, upstream = 200, downstream = 200, BSgenomeName)
{
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
    if ( ! file.exists(outputDir))
    {
        dir.create(outputDir)
    }
    OfftargetFile <-paste(outputDir, "OfftargetAnalysis.xls", sep = "")
    OfftargetSummary <-paste(outputDir, "Summary.xls", sep = "")
    gRNAsPlusPAM<- unique(scores$gRNAPlusPAM)
    names <- gRNAsPlusPAM
    top5OfftargetTotalScore <- numeric(length(names))
    topNOfftargetTotalScore <- top5OfftargetTotalScore
    temp <- cbind(names, gRNAsPlusPAM, top5OfftargetTotalScore, 
        topNOfftargetTotalScore)
    mismatche.distance2PAM <- matrix(ncol = 11, nrow = length(names))
    append <- FALSE
    for (i in 1:length(gRNAsPlusPAM))
    {
        this.score <- scores[scores$gRNAPlusPAM == gRNAsPlusPAM[i],]
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
        mismatche.distance2PAM[i,] <- 
            as.character(forSummary$mismatche.distance2PAM)
	forSummary <- forSummary[2:11,]
	forSummary <- forSummary[order(forSummary$score, decreasing=TRUE),]
	mismatche.distance2PAM[i,2:11] <-
            as.character(forSummary$mismatche.distance2PAM)
        if (annotateExon)
        {
            	score.RD <- GRanges(seqnames = Rle(this.score$chrom), 
                	ranges = IRanges(start = this.score$chromStart, 
                	end = this.score$chromEnd, names = this.score$forViewInUCSC))
            	allExons <- as(exons(txdb, columns="gene_id"),"GRanges")
            	allExons <- allExons[as.character(seqnames(allExons)) %in% 
                	unique(as.character(seqnames(score.RD))),]
            	ann.scores <- overlapsAny(score.RD, allExons, minoverlap = 1L, 
                	type = "any",ignore.strand=TRUE)
		introns =  unlist(intronsByTranscript(txdb))
		introns <- introns[as.character(seqnames(introns)) %in% 
		unique(as.character(seqnames(score.RD))),]
		ann.scores.intron <- overlapsAny(score.RD, introns, minoverlap = 1L, 
			type = "any", ignore.strand=TRUE)
		inIntron <- cbind(forViewInUCSC = names(score.RD),
			inIntron = unlist(ann.scores.intron))
        	inIntron[inIntron[,2] == FALSE, 2] <- ""
            	inExon <- cbind(forViewInUCSC = names(score.RD),
                	inExon = unlist(ann.scores))
            	inExon[inExon[,2] == FALSE, 2] <- ""
		overlapGenes <- findOverlaps(score.RD, allExons, minoverlap = 1L, 
			type = "any",ignore.strand=TRUE)
		overlapGenes.ID <- unlist(allExons[subjectHits(overlapGenes),]$gene_id)
		this.score <- cbind(this.score, entrez_id = "", symbol="")
		overlapGenes.ID <- subset(overlapGenes.ID, !is.na(overlapGenes.ID))
		query.ind <- subset(queryHits(overlapGenes), !is.na(overlapGenes.ID))
		if (length(overlapGenes.ID) > 0 && !missing(orgAnn) )
		{
			print(overlapGenes.ID)
			overlapGenes.symbol <- as.vector(unlist(mget(overlapGenes.ID, 
				envir= orgAnn, ifnotfound=NA)))
			this.score$symbol = as.character(this.score$symbol)
			this.score$symbol[query.ind] = overlapGenes.symbol
			this.score$symbol[is.na(this.score$entrez_id)] = ""
		}
		this.score$entrez_id = as.character(this.score$entrez_id)
		this.score$entrez_id[query.ind] = overlapGenes.ID
        	this.score <- merge(this.score,inExon)
		this.score <- merge(this.score,inIntron)
		this.score$entrez_id[is.na(this.score$entrez_id)] = ""
		if (length(grep("FBgn",overlapGenes.ID[1])) > 0 && !missing(orgAnn))
			{
				
				temp.id <- this.score$entrez_id
				this.score$entrez_id <- this.score$symbol
				this.score$symbol <- temp.id
				rm(temp.id)
			}
		if (!missing(orgAnn))
		{
                  this.score <- cbind(name = this.score$name,
               		gRNAPlusPAM = this.score$gRNAPlusPAM,
                	OffTargetSequence = this.score$OffTargetSequence,
                	inExon = as.character(this.score$inExon),
			inIntron = as.character(this.score$inIntron),
			entrez_id = this.score$entrez_id,
			gene = this.score$symbol,
                	score = this.score$score, n.mismatch = this.score$n.mismatch, 
                	mismatche.distance2PAM = as.character(
                	this.score$mismatche.distance2PAM), 
                	alignment = this.score$alignment,
                	NGG = as.character(this.score$NGG),
                	toViewInUCSC = this.score$forViewInUCSC, 
                	strand = this.score$strand,
                	chrom = this.score$chrom, chromStart = this.score$chromStart,
                	chromEnd = this.score$chromEnd)
		}
		else
		{
		    this.score <- cbind(name = this.score$name,
                        gRNAPlusPAM = this.score$gRNAPlusPAM,
                        OffTargetSequence = this.score$OffTargetSequence,
                        inExon = as.character(this.score$inExon),
                        inIntron = as.character(this.score$inIntron),
                        entrez_id = this.score$entrez_id,
                        score = this.score$score, n.mismatch = this.score$n.mismatch,
                        mismatche.distance2PAM = as.character(
                        this.score$mismatche.distance2PAM),
                        alignment = this.score$alignment,
                        NGG = as.character(this.score$NGG),
                        toViewInUCSC = this.score$forViewInUCSC,
                        strand = this.score$strand,
                        chrom = this.score$chrom, chromStart = this.score$chromStart,
                        chromEnd = this.score$chromEnd)
		}
        }
        else
            this.score <- cbind(name = this.score$name, 
                gRNAPlusPAM = this.score$gRNAPlusPAM,
                OffTargetSequence = this.score$OffTargetSequence,
                score = this.score$score, n.mismatch = this.score$n.mismatch, 
                mismatche.distance2PAM = 
                as.character(this.score$mismatche.distance2PAM), 
                alignment = this.score$alignment,
                NGG = as.character(this.score$NGG),
                toViewInUCSC = this.score$forViewInUCSC, 
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
    temp <- cbind(temp, mismatche.distance2PAM)
    colnames(temp)[5] <- "top1Hit(onTarget)MMdistance2PAM"
    #colnames(temp)[3] <- paste("top", min(5, maxN - 1), 
    #    "OfftargetTotalScore", sep = "")
    #colnames(temp)[4] <- paste("top", maxN.totalScore - 1, 
    #    "OfftargetTotalScore", sep = "")
    colnames(temp)[4] <- paste("top", topN.OfftargetTotalScore, 
       "OfftargetTotalScore", sep = "")
    colnames(temp)[6:15] <- paste("topOfftarget", 1:10, "MMdistance2PAM",
        sep = "")
    #if (maxN == maxN.totalScore)
    #    temp <- temp[,-4]
    Offtargets <- read.table(OfftargetFile, sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE)
    if (fetchSequence)
    {
        Start <- as.numeric(Offtargets$chromStart) - as.numeric(upstream)
        End <- as.numeric(Offtargets$chromEnd) + as.numeric(downstream)
        strand <- Offtargets$strand		
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
    if (annotateExon)
    {
    	Offtargets$inExon[is.na(Offtargets$inExon)] <- ""
    	Offtargets$inIntron[is.na(Offtargets$inIntron)] <- ""
    	Offtargets$entrez_id[is.na(Offtargets$entrez_id)] <- ""
    }
    write.table(Offtargets[order(as.character(Offtargets$name), 
        -as.numeric(as.character(Offtargets$score)), 
        as.character(Offtargets$OffTargetSequence)),], 
        file = OfftargetFile, sep = "\t", row.names = FALSE)
    list(offtargets = Offtargets,summary = temp)
}
