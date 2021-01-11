filterOffTargetWithoutBSgenome <- function (scores, 
    min.score = 0.01, topN = 200, topN.OfftargetTotalScore = 20, 
    annotateExon = TRUE, txdb, orgAnn, ignore.strand = TRUE, 
    outputDir, oneFilePergRNA = FALSE, fetchSequence = TRUE, 
    upstream = 200, downstream = 200, BSgenomeName, baseBeforegRNA = 4, 
    baseAfterPAM = 3, featureWeightMatrixFile = system.file("extdata", 
        "DoenchNBT2014.csv", package = "CRISPRseek"), rule.set = c("Root_RuleSet1_2014", 
        "Root_RuleSet2_2016", "CRISPRscan"),
    genomeSeq, useBSgenome = FALSE) 
{
    rule.set <- match.arg(rule.set)
    if (featureWeightMatrixFile != "" && file.exists(featureWeightMatrixFile)) {
        featureWeightMatrix <- read.csv(featureWeightMatrixFile, 
            header = TRUE)
    }
    if (fetchSequence && (missing(BSgenomeName) || class(BSgenomeName) != 
        "BSgenome") && useBSgenome) {
        stop("To fetch sequences, BSgenomeName is required as BSgenome object!")
    }
    if (!useBSgenome && (missing(genomeSeq) || class(genomeSeq) != "DNAStringSet"))
        stop("To fetch sequence from non-BSgenome, genomeSeq needs to be DNAStringSet class!")
    if (annotateExon && (missing(txdb) || (class(txdb) != "TxDb" && 
        class(txdb) != "TranscriptDb"))) {
        stop("To indicate whether an offtarget is inside an exon, txdb is\n required as TxDb object!")
    }
    if (annotateExon && missing(orgAnn)) {
        warning("orgAnn was not included. See the updated manual for information about how to use the orgAnn parameter to generate gene identifiers in the offTarget output file")
    }
    scores <- scores[scores$score >= min.score, ]
    if (length(grep("IsMismatch.pos", colnames(scores))) > 0) 
        scores <- scores[, -c(grep("IsMismatch.pos", colnames(scores)))]
    if (substr(outputDir, nchar(outputDir), nchar(outputDir)) != 
        .Platform$file.sep) {
        outputDir <- paste(outputDir, "", sep = .Platform$file.sep)
    }
    if (!file.exists(outputDir)) {
        dir.create(outputDir)
    }
    OfftargetFile <- paste(outputDir, "OfftargetAnalysis.xls", 
        sep = "")
    OfftargetSummary <- paste(outputDir, "Summary.xls", sep = "")
    gRNAsPlusPAM <- unique(scores$name)
    names <- gRNAsPlusPAM
    top5OfftargetTotalScore <- numeric(length(names))
    topNOfftargetTotalScore <- top5OfftargetTotalScore
    temp <- cbind(names, gRNAsPlusPAM, top5OfftargetTotalScore, 
        topNOfftargetTotalScore)
    mismatch.distance2PAM <- matrix(ncol = 11, nrow = length(names))
    append <- FALSE
    for (i in 1:length(gRNAsPlusPAM)) {
        this.score <- scores[scores$name == gRNAsPlusPAM[i], 
            ]
        this.score <- this.score[order(this.score$score, this.score$n.mismatch, 
            decreasing = c(TRUE, FALSE)), ]
        maxN <- min(topN + 1, dim(this.score)[1])
        this.score <- this.score[1:maxN, ]
        maxN.totalScore <- min(maxN, (topN.OfftargetTotalScore + 
            1))
        if (this.score$n.mismatch[1] == 0 && as.numeric(as.character(this.score$NGG[1])) == 
            1) {
            start.ind <- 2
            end.ind <- min(maxN, 6)
            end.forSummary <- 11
        }
        else {
            start.ind <- 1
            maxN <- maxN - 1
            maxN.totalScore <- maxN.totalScore - 1
            end.forSummary <- 10
            end.ind <- min(maxN, 5)
        }
        temp[i, 3] <- sum(this.score$score[start.ind:end.ind])
        if (maxN < maxN.totalScore) 
            temp[i, 4] <- sum(this.score$score[start.ind:maxN])
        else temp[i, 4] <- sum(this.score$score[start.ind:maxN.totalScore])
        temp[i, 2] <- unique(this.score$gRNAPlusPAM)
        mismatch.distance2PAM[i, ] <- ifelse(as.character(this.score$mismatch.distance2PAM[1]) == 
            "", "NMM", "perfect match not found")
        forSummary <- this.score[start.ind:end.forSummary, ]
        forSummary <- forSummary[order(forSummary$score, decreasing = TRUE), 
            ]
        mismatch.distance2PAM[i, 2:11] <- as.character(forSummary$mismatch.distance2PAM)
        if (dim(forSummary)[1] < 10) 
            mismatch.distance2PAM[i, (dim(forSummary)[1] + 1):11] <- "NA"
        this.score <- cbind(name = this.score$name, gRNAPlusPAM = this.score$gRNAPlusPAM, 
            OffTargetSequence = this.score$OffTargetSequence, 
            score = this.score$score, n.mismatch = this.score$n.mismatch, 
            mismatch.distance2PAM = as.character(this.score$mismatch.distance2PAM), 
            alignment = this.score$alignment, NGG = as.character(this.score$NGG), 
            forViewInUCSC = this.score$forViewInUCSC, strand = this.score$strand, 
            chrom = this.score$chrom, chromStart = this.score$chromStart, 
            chromEnd = this.score$chromEnd)
        if (oneFilePergRNA & dim(this.score)[1] > 0) 
            write.table(this.score[!is.na(this.score[, grep("score", 
                colnames(this.score))]), ], file = paste(outputDir, 
                "OfftargetAnalysis-", as.character(temp[i, 1]), 
                ".xls", sep = ""), sep = "\t", row.names = FALSE)
        if (i == 1 && dim(this.score)[1] > 0) {
            write.table(this.score, file = OfftargetFile, sep = "\t", 
                row.names = FALSE, append = append)
            append <- TRUE
        }
        else if (dim(this.score)[1] > 0) {
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
    if (annotateExon) {
        Offtargets <- annotateOffTargets(Offtargets, txdb, orgAnn, 
            ignore.strand)
    }
    ontargets <- subset(Offtargets, Offtargets$n.mismatch == 
        0)
    if (dim(ontargets)[1] > 0) {
        chr <- as.character(ontargets$chrom)
        strand <- as.character(ontargets$strand)
        Start <- ifelse(strand == "-", as.numeric(as.character(ontargets$chromStart)) - 
            baseAfterPAM, as.numeric(as.character(ontargets$chromStart)) - 
            baseBeforegRNA)
        End <- ifelse(strand == "-", as.numeric(as.character(ontargets$chromEnd)) + 
            as.numeric(baseBeforegRNA), as.numeric(as.character(ontargets$chromEnd)) + 
            as.numeric(baseAfterPAM))
        starts <- unlist(apply(cbind(Start, 1), 1, max))
        if (useBSgenome)
        {
          ends <- unlist(apply(cbind(End, seqlengths(BSgenomeName)[chr]),
            1, min))
          extendedSequence <- getSeq(BSgenomeName, names = chr, 
            start = starts, end = ends, strand = strand, width = NA, 
            as.character = TRUE)
        }
        else
        {
            ends <- unlist(apply(cbind(End, width(genomeSeq)[names(genomeSeq) %in% chr]),
                1, min))
            extended.info <- data.frame(chrom = chr, start = starts, end = ends, strand = strand)
            extendedSequence <- getSeq(genomeSeq, as(extended.info, "GRanges"))
        }
          
        ontargets <- cbind(ontargets, extendedSequence = extendedSequence)
        if (rule.set == "Root_RuleSet1_2014") {
            gRNAefficiency <- calculategRNAEfficiency(extendedSequence, 
                baseBeforegRNA = baseBeforegRNA, featureWeightMatrix = featureWeightMatrix)
        }
        else if (rule.set == "Root_RuleSet2_2016") {
            gRNAefficiency <- calculategRNAEfficiency2(extendedSequence)
        }
        else if (rule.set == "CRISPRscan") {
            gRNAefficiency <- calculategRNAEfficiencyCRISPRscan(extendedSequence, 
                featureWeightMatrix = featureWeightMatrix)
        }
        ontargets <- cbind(ontargets, gRNAefficacy = gRNAefficiency)
        Offtargets <- merge(Offtargets, ontargets, all = TRUE)
    }
    if (fetchSequence) {
        strand <- as.character(Offtargets$strand)
        chr <- as.character(Offtargets$chrom)
        Start <- ifelse(strand == "-", as.numeric(as.character(Offtargets$chromStart)) - 
            as.numeric(downstream), as.numeric(as.character(Offtargets$chromStart)) - 
            as.numeric(upstream))
        End <- ifelse(strand == "-", as.numeric(as.character(Offtargets$chromEnd)) + 
            as.numeric(upstream), as.numeric(as.character(Offtargets$chromEnd)) + 
            as.numeric(downstream))
        starts <- unlist(apply(cbind(Start, 1), 1, max))
        if (useBSgenome)
        {
            ends <- unlist(apply(cbind(End, seqlengths(BSgenomeName)[chr]),
                1, min))
            seq <- getSeq(BSgenomeName, names = chr, start = starts, 
               end = ends, strand = strand, width = NA, as.character = TRUE)
        }
        else
        {
            #ends <- unlist(apply(cbind(End, width(genomeSeq)[names(genomeSeq) %in% chr]),
            chrLen <- unlist(lapply(chr, function(i) width(genomeSeq)[names(genomeSeq) == i]))
            ends <- unlist(apply(cbind(End, chrLen),
                1, min))
            extended.info <- data.frame(chrom = chr, start = starts, end = ends, strand = strand)            
       	    seq <- getSeq(genomeSeq, as(extended.info, "GRanges")) 
        }
        Offtargets <- cbind(Offtargets, flankSequence = seq)
    }
    colnames(Offtargets)[colnames(Offtargets) == "NGG"] = "isCanonicalPAM"
    write.table(temp, file = OfftargetSummary, sep = "\t", row.names = FALSE)
    write.table(Offtargets[order(as.character(Offtargets$name), 
        -as.numeric(as.character(Offtargets$score)), as.character(Offtargets$OffTargetSequence)), 
        ], file = OfftargetFile, sep = "\t", row.names = FALSE)
    list(offtargets = unique(Offtargets), summary = unique(temp))
}

