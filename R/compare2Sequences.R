compare2Sequences <- function(inputFile1Path, inputFile2Path, inputNames=c("Seq1", "Seq2"), 
	format = c("fasta", "fasta"), header = FALSE, findgRNAsWithREcutOnly = FALSE,
    searchDirection = c("both","1to2", "2to1"), BSgenomeName,
    REpatternFile=system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek"),
    minREpatternSize = 6, findgRNAs = c(TRUE, TRUE), removegRNADetails = c(FALSE, FALSE), 
    exportAllgRNAs = c("no", "all", "fasta", "genbank"), annotatePaired =  FALSE,
    overlap.gRNA.positions = c(17, 18), findPairedgRNAOnly = FALSE, 
    min.gap = 0, max.gap = 20, gRNA.name.prefix = "_gR", PAM.size = 3, 
    gRNA.size = 20, PAM = "NGG", PAM.pattern = "N[A|G]G$",
    allowed.mismatch.PAM = 1, max.mismatch = 3, 
    outputDir, upstream =0, downstream = 0,
    weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 
    0.685, 0.583), overwrite = FALSE, baseBeforegRNA = 4, 
    baseAfterPAM = 3, featureWeightMatrixFile = system.file("extdata", 
       "DoenchNBT2014.csv", package = "CRISPRseek"), foldgRNAs = FALSE, 
        gRNA.backbone="GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU",
        temperature = 37,
        scoring.method = c("Hsu-Zhang", "CFDscore"),
        subPAM.activity = hash( AA =0,
          AC =   0,
          AG = 0.259259259,
          AT = 0,
          CA = 0,
          CC = 0,
          CG = 0.107142857,
          CT = 0,
          GA = 0.069444444,
          GC = 0.022222222,
          GG = 1,
          GT = 0.016129032,
          TA = 0,
          TC = 0,
          TG = 0.038961039,
          TT = 0),
     subPAM.position = c(22, 23),
     PAM.location = "3prime",
     mismatch.activity.file = system.file("extdata", 
         "NatureBiot2016SuppTable19DoenchRoot.csv", 
         package = "CRISPRseek")
    )
{
        scoring.method <- match.arg(scoring.method)
        exportAllgRNAs <- match.arg(exportAllgRNAs)
        searchDirection <- match.arg(searchDirection)
        if (scoring.method == "Hsu-Zhang")
        {
             if (length(weights) !=  gRNA.size)
                 stop("Please make sure the size of weights vector 
                     equals to the gRNA.size!\n")
        }
        else if (scoring.method ==  "CFDscore")
        {
            mismatch.activity <- read.csv(mismatch.activity.file)
            required.col <- c("Mismatch.Type", "Position", "Percent.Active")
            if (length(intersect(colnames(mismatch.activity), required.col)) !=
                length(required.col))
                stop("Please rename the mismatch activity file column to contain at least
                   these 3 column names: Mismatch.Type, Position, Percent.Active\n")
        }
        if (class(inputFile1Path) != "DNAStringSet" || class(inputFile2Path) != "DNAStringSet")
        {
	    if ((format[1] == "bed" || format[2] == "bed") && 
                (missing(BSgenomeName) || class(BSgenomeName) != "BSgenome"))
                   stop("BSgenomeName is required as BSgenome object when input file is in bed format!")
        }
	append = ifelse(overwrite, FALSE, TRUE)
	if (class(inputFile1Path) != "DNAStringSet")
	{
	    outputDir1 <- file.path(outputDir, paste(basename(inputFile1Path),
               format(Sys.time(), "%b-%d-%Y"), sep="-"))
	}
	else
	{
           outputDir1 <- file.path(outputDir, "File1")
	}
	if (class(inputFile2Path) != "DNAStringSet")
	{
	   outputDir2 <- file.path(outputDir, paste(basename(inputFile2Path),
               format(Sys.time(), "%b-%d-%Y"), sep="-"))
	}
	else
	{
	   outputDir2 <- file.path(outputDir, "File2")
	}
	if(searchDirection == "both" || searchDirection == "1to2")
	{
		cat("search for gRNAs for input file1...\n")
                if (findgRNAs[1])
                {
		    tryCatch(
			 (gRNAs1 = offTargetAnalysis(inputFile1Path, format = format[1], 
			     findgRNAs = findgRNAs[1], annotatePaired =  annotatePaired,
                             exportAllgRNAs = exportAllgRNAs, gRNAoutputName= inputNames[1],
			     findPairedgRNAOnly = findPairedgRNAOnly, chromToSearch = "", 
			     findgRNAsWithREcutOnly = findgRNAsWithREcutOnly, 
			     REpatternFile = REpatternFile, minREpatternSize = minREpatternSize, 
			     overlap.gRNA.positions =  overlap.gRNA.positions,
			     min.gap = min.gap, max.gap = max.gap, 
			     gRNA.name.prefix = gRNA.name.prefix, PAM.size = PAM.size,
			     gRNA.size = gRNA.size, PAM = PAM, PAM.pattern = PAM.pattern,
                             outputUniqueREs = FALSE,
			     outputDir = outputDir1, upstream.search = upstream,
			     downstream.search = downstream,
			     weights = weights, foldgRNAs = FALSE, overwrite = overwrite,
			     featureWeightMatrixFile = featureWeightMatrixFile, 
            		     baseBeforegRNA = baseBeforegRNA, BSgenomeName = BSgenomeName,
            		     baseAfterPAM = baseAfterPAM, header = header,
                             subPAM.position = subPAM.position,
                             subPAM.activity = subPAM.activity,
                             PAM.location = PAM.location,
                             mismatch.activity.file = mismatch.activity.file)), 
			 error = function(e) {print(e); gRNAs1 = DNAStringSet()})
             }
             else if (class(inputFile1Path) == "DNAStringSet") 
             {
                gRNAs1 <- inputFile1Path
             }
            else  if (format[1] == "fasta" || format[1] == "fastq")
            {
                gRNAs1 <- readDNAStringSet(inputFile1Path, format[1],
                      use.names = TRUE)
            }
            else
            {
                stop("format needs to be either fasta or fastq for gRNA file!")
            }
	}
	if(searchDirection == "both" || searchDirection == "2to1")
	{
             if(findgRNAs[2])
             {
		cat("search for gRNAs for input file2...\n")
		tryCatch((gRNAs2 = offTargetAnalysis(inputFile2Path, format = format[2],          
			findgRNAs = findgRNAs[2], annotatePaired =  annotatePaired,
                        exportAllgRNAs = exportAllgRNAs, gRNAoutputName = inputNames[2],
			findPairedgRNAOnly = findPairedgRNAOnly, chromToSearch = "",
            		findgRNAsWithREcutOnly = findgRNAsWithREcutOnly, 
            		REpatternFile = REpatternFile, minREpatternSize = minREpatternSize,
                        outputUniqueREs = FALSE,
            		overlap.gRNA.positions =  overlap.gRNA.positions, 
            		min.gap = min.gap, max.gap = max.gap, 
            		gRNA.name.prefix = gRNA.name.prefix, PAM.size = PAM.size,
            		gRNA.size = gRNA.size, PAM = PAM, PAM.pattern = PAM.pattern, 
            		outputDir = outputDir2,upstream.search = upstream,
			downstream.search = downstream, 
            		weights = weights, foldgRNAs = FALSE, overwrite = overwrite,
                        featureWeightMatrixFile = featureWeightMatrixFile,
                        baseBeforegRNA = baseBeforegRNA, BSgenomeName = BSgenomeName,
                        baseAfterPAM = baseAfterPAM, header = header, 
                        subPAM.position = subPAM.position,
                        subPAM.activity = subPAM.activity,
                        PAM.location = PAM.location,
                        mismatch.activity.file = mismatch.activity.file)), 
			error=function(e) {print(e); gRNAs2 = DNAStringSet()})
            }
            else if (class(inputFile2Path) == "DNAStringSet")
            {
                gRNAs2 <- inputFile2Path
            }
            else  if (format[2] == "fasta" || format[2] == "fastq")
            {
                gRNAs2 <- readDNAStringSet(inputFile2Path, format[2],
                      use.names = TRUE)
            }
            else
            {
                stop("format needs to be either fasta or fastq for gRNA file!")
            }
	}
    print("Scoring ...")
    if (class(inputFile1Path) != "DNAStringSet")
    {
	if (format[1] == "bed")
	    subjects1 <- getSeqFromBed(inputFile1Path, header = header,
                BSgenomeName = BSgenomeName, upstream = upstream, downstream = downstream)
	else
    	    subjects1 <- readDNAStringSet(inputFile1Path, format = format[1],
    	        use.names = TRUE)
    }
    else
    {
	   subjects1 <- inputFile1Path
    }
   if (class(inputFile2Path) != "DNAStringSet")
   {
	if (format[2] == "bed")
            subjects2 <- getSeqFromBed(inputFile2Path, header = header, 
	        BSgenomeName = BSgenomeName, upstream = upstream, downstream = downstream)
	else
    	    subjects2 <- readDNAStringSet(inputFile2Path, format = format[2],
       	        use.names = TRUE)
   }
   else
   {
	   subjects2 <- inputFile2Path
   }
    outfile <- tempfile(tmpdir = getwd())
    seqname <- names(subjects2)
    seqname <- gsub("'", "", seqname)
    seqname <- gsub(" ", "", seqname)
    seqname <- gsub("\t", ":", seqname)
    names(subjects2) <- seqname
    #write.table(subjects2, file="testSeqs")
    #revsubject <- reverseComplement(subjects2)
    #chrom.len <- nchar(as.character(subjects2))
    if(searchDirection == "both" || searchDirection == "1to2")
	{
	   for (j in 1:length(subjects2))
	   {
                if (j == 1)
                    hits <- searchHits(gRNAs = gRNAs1, PAM = PAM, PAM.pattern = PAM.pattern, 
                        seqs = subjects2[[j]], seqname = names(subjects2)[j],
                        max.mismatch = max.mismatch, PAM.size = PAM.size, 
                        gRNA.size = gRNA.size, allowed.mismatch.PAM = allowed.mismatch.PAM,
                        PAM.location = PAM.location, outfile = outfile) 
                else
                    hits <- rbind(hits, searchHits(gRNAs = gRNAs1, PAM = PAM, 
                        PAM.pattern = PAM.pattern,
                        seqs = subjects2[[j]], seqname = names(subjects2)[j],
                        max.mismatch = max.mismatch, PAM.size = PAM.size,
                        gRNA.size = gRNA.size, allowed.mismatch.PAM = allowed.mismatch.PAM,
                        PAM.location = PAM.location, outfile = outfile))
           }
	} # end of if searchDirection == "both" or searchDirection == "1to2"
	cat("finish off-target search in sequence 2\n") 
   	seqname <- names(subjects1)
	seqname <- gsub("'", "", seqname)
	seqname <- gsub(" ", "", seqname)
	seqname <- gsub("\t", ":", seqname)
	#revsubject <- reverseComplement(subjects1)
	#revsubject <- reverseComplement(subjects1[[1]])
	names(subjects1) <- seqname
	#chrom.len = nchar(as.character(subjects1))
	if(searchDirection == "both" || searchDirection == "2to1")
	{
	    for (j in 1:length(subjects1))
	    {
                if (j == 1  && !exists("hits"))
                    hits <- searchHits(gRNAs = gRNAs2, PAM = PAM, PAM.pattern = PAM.pattern, 
                       seqs = subjects1[[j]], seqname = names(subjects1)[j],
                       max.mismatch = max.mismatch, PAM.size = PAM.size, 
                       gRNA.size = gRNA.size, allowed.mismatch.PAM = allowed.mismatch.PAM,
                       PAM.location = PAM.location, outfile = outfile) 
               else
                   hits <- rbind(hits, searchHits(gRNAs = gRNAs2, PAM = PAM, 
                       PAM.pattern = PAM.pattern,
                       seqs = subjects1[[j]], seqname = names(subjects1)[j],
                       max.mismatch = max.mismatch, PAM.size = PAM.size, 
                       gRNA.size = gRNA.size, allowed.mismatch.PAM = allowed.mismatch.PAM,
                       PAM.location = PAM.location, outfile = outfile))
	    }
	} # if searchDirection == "both" or searchDirection == "2to1"
	cat("finish off-target search in sequence 1\n")
        if (removegRNADetails[1])
            unlink(outputDir1, recursive = TRUE)
        if (removegRNADetails[2])
            unlink(outputDir2, recursive = TRUE)
	featureVectors <- buildFeatureVectorForScoring(hits = hits, 
	    canonical.PAM = PAM, gRNA.size = gRNA.size,
            subPAM.position = subPAM.position, PAM.location = PAM.location)
	cat("finish feature vector building\n")
        if ( scoring.method ==  "CFDscore")
            scores <- getOfftargetScore2(featureVectors,
                subPAM.activity = subPAM.activity,
                mismatch.activity.file = mismatch.activity.file)
        else
	    scores <- getOfftargetScore(featureVectors, weights = weights)
	cat("finish score calculation\n")
	targetInSeq1 <- scores$gRNAPlusPAM
	targetInSeq2 <- scores$gRNAPlusPAM
        if ( scoring.method ==  "CFDscore")
        {	
            scoreForSeq1 <- rep(1, dim(scores)[1])
	    scoreForSeq2 <- rep(1, dim(scores)[1])
            max.score <- 1
        }
        else
        {
            scoreForSeq1 <- rep(100, dim(scores)[1])
            scoreForSeq2 <- rep(100, dim(scores)[1])
            max.score <- 100
        }
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
	else if (searchDirection == "1to2")
	{
		scoreForSeq2 = scores$score
		targetInSeq2 = scores$OffTargetSequence
		targetSeqName = scores$name
	}
	else if (searchDirection == "2to1")
	{
		scoreForSeq1 = scores$score
		targetInSeq1 = scores$OffTargetSequence
		targetSeqName = scores$name
	}
	else
	{	
		seqNames <- c(names(subjects1), names(subjects2))
		fileIndex <- c(rep(1, length(subjects1)), rep(2, length(subjects2)))
		offTargetFiles <- fileIndex[match(scores$chrom,seqNames)]
		targetSeqName <- scores$name
		scoreForSeq1[offTargetFiles == 1] <- 
			scores[offTargetFiles == 1, ]$score
		scoreForSeq2[offTargetFiles == 2] <- 
			scores[offTargetFiles == 2, ]$score
		targetInSeq1[offTargetFiles ==1] <- 
			scores[offTargetFiles == 1,]$OffTargetSequence
		targetInSeq2[offTargetFiles == 2] <- 
			scores[offTargetFiles == 2,]$OffTargetSequence
	}
	targetSeqName <- gsub(paste(gRNA.name.prefix, "[0-9]+[f|r]", sep=""), "", targetSeqName)
	targetSeqName <- gsub(paste("_gR", "[0-9]+[f|r]", sep=""), "", targetSeqName)
	targetSeqName <- gsub( "Start[0-9]+End[0-9]+", "", targetSeqName)		
	seqs.new <- as.data.frame(cbind(name = scores$name,
		gRNAPlusPAM = scores$gRNAPlusPAM,
		targetInSeq1 = targetInSeq1,
		targetInSeq2 = targetInSeq2,
	    guideAlignment2OffTarget = scores$alignment,
	    offTargetStrand = scores$strand,
		scoreForSeq1 = scoreForSeq1,
		scoreForSeq2 = scoreForSeq2,
		mismatch.distance2PAM = as.character(scores$mismatch.distance2PAM),
		n.mismatch = scores$n.mismatch,
		offTarget = scores$forViewInUCSC,
		targetSeqName = targetSeqName)
	)
	 
	if(searchDirection == "both" || searchDirection == "1to2")
	{
	if (length(setdiff(names(gRNAs1),seqs.new[,1])) >0)
	{
		gRNAnames <- names(gRNAs1)[!names(gRNAs1) %in% seqs.new[,1]]
		targetSeqName <- gsub(paste(gRNA.name.prefix, "[0-9]+[f|r]", sep=""), "", gRNAnames)
		targetSeqName <- gsub(paste("_gR", "[0-9]+[f|r]", sep=""), "", targetSeqName)
		targetSeqName <- gsub( "Start[0-9]+End[0-9]+", "", targetSeqName)		
		seqs1.only <- as.data.frame(cbind(name = gRNAnames,
			gRNAPlusPAM = as.character(gRNAs1)[names(gRNAs1) %in% gRNAnames],
			targetInSeq1 = as.character(gRNAs1)[names(gRNAs1) %in% gRNAnames],
			targetInSeq2 = rep("NA", length(targetSeqName)),
		  guideAlignment2OffTarget = rep("NA", length(targetSeqName)),
		  offTargetStrand = rep("NA", length(targetSeqName)),
			scoreForSeq1 = rep(max.score, length(targetSeqName)),
			scoreForSeq2 = rep(0, length(targetSeqName)),
			mismatch.distance2PAM = rep("NA", length(targetSeqName)),
			n.mismatch = rep("NA", length(targetSeqName)),
			offTarget = rep("NA", length(targetSeqName)),
			targetSeqName = targetSeqName)
			)
		seqs.new <- rbind(seqs.new, seqs1.only)
	}
	}
	if(searchDirection == "both" || searchDirection == "2to1")
	{
	if (length(setdiff(names(gRNAs2),seqs.new[,1])) >0)
	{
		gRNAnames <- names(gRNAs2)[!names(gRNAs2) %in% seqs.new[,1]]
		targetSeqName <- gsub(paste(gRNA.name.prefix, "[0-9]+[f|r]", sep=""), "",gRNAnames)
		targetSeqName <- gsub(paste("_gR", "[0-9]+[f|r]", sep=""), "", targetSeqName)
		targetSeqName <- gsub( "_Start[0-9]+End[0-9]+", "", targetSeqName)
		seqs2.only <- as.data.frame(cbind(name = gRNAnames,
    		gRNAPlusPAM = as.character(gRNAs2)[names(gRNAs2) %in% gRNAnames],
		    targetInSeq1 = rep("NA", length(targetSeqName)),
			  targetInSeq2 = as.character(gRNAs2)[names(gRNAs2) %in% gRNAnames],
		    guideAlignment2OffTarget = rep("NA", length(targetSeqName)),
        offTargetStrand = rep("NA", length(targetSeqName)),
        scoreForSeq1 = rep(0, length(targetSeqName)),
			  scoreForSeq2 = rep(max.score, length(targetSeqName)),
			  mismatch.distance2PAM = rep("NA", length(targetSeqName)),
			  n.mismatch = rep("NA", length(targetSeqName)),
			  offTarget = rep("NA", length(targetSeqName)),
			  targetSeqName = targetSeqName)
		)
		seqs.new <- rbind(seqs.new, seqs2.only)
	}
	}
	seqs <- unique(cbind(seqs.new, gRNAefficacy = 0, 
        scoreDiff = round(as.numeric(as.character(seqs.new[,7])) - 
                            as.numeric(as.character(seqs.new[,8]),4), 3)))
#	rownames(seqs) <- seqs[,1]
#	seqs <- as.data.frame(seqs)
#save(seqs, file="seqs")
        if (substr(outputDir1, nchar(outputDir1), nchar(outputDir1)) != .Platform$file.sep)
    	{
       	    outputDir1 <- paste(outputDir1, "", sep = .Platform$file.sep)
    	}
	if (substr(outputDir2, nchar(outputDir2), nchar(outputDir2)) != .Platform$file.sep)
        {
            outputDir2 <- paste(outputDir2, "", sep = .Platform$file.sep)
        }
	if (searchDirection == "both")
	{
            if (findgRNAs[1])
            {
		eff1File <- paste(outputDir1, "gRNAefficacy.xls", sep = "")
                gRNAeff1 <- read.table(eff1File,sep="\t", header=TRUE, stringsAsFactors=FALSE)
            }  
            if (findgRNAs[2])
            {
		eff2File <- paste(outputDir2, "gRNAefficacy.xls", sep = "")
		gRNAeff2 <- read.table(eff2File,sep="\t", header=TRUE, stringsAsFactors=FALSE)
            }
            if (findgRNAs[1])
            {
                if (findgRNAs[2])
		    gRNAeff <- rbind(gRNAeff1, gRNAeff2) 
                else
                    gRNAeff <- gRNAeff1
            }
            else if (findgRNAs[2])
                gRNAeff <- gRNAeff2
        }
	if(searchDirection == "1to2" && findgRNAs[1])
	{
		eff1File <- paste(outputDir1, "gRNAefficacy.xls", sep = "")
		gRNAeff <- read.table(eff1File,sep="\t", header=TRUE, stringsAsFactors=FALSE)
        }
	if(searchDirection == "2to1" && findgRNAs[2])
	{
		eff2File <- paste(outputDir2, "gRNAefficacy.xls", sep = "")
		gRNAeff <- read.table(eff2File,sep="\t", header=TRUE, stringsAsFactors=FALSE)
	}
        if(exists("gRNAeff"))
        {
	    m <- match(seqs$name, gRNAeff$name)
	    seqs$gRNAefficacy <- gRNAeff$gRNAefficacy[m]
        }
	originalDir <- getwd()
	setwd(outputDir)
        if (foldgRNAs)
        {
           gRNAs.withoutPAM <- substr(as.character(seqs$gRNAPlusPAM), 1, gRNA.size)
           folded.gRNAs <- foldgRNAs(gRNAs.withoutPAM, gRNA.backbone = gRNA.backbone, 
           temperature = temperature)
           if (length(dim(folded.gRNAs)) > 0)
               seqs <- cbind(seqs, folded.gRNAs[,-1])
        }
	if (dim(seqs)[1] ==1)
	{
		write.table(seqs, file = "scoresFor2InputSequences.xls",
			sep = "\t", row.names = FALSE, col.names=TRUE)
	}
	else
	{
             write.table(seqs[order(as.numeric(seqs[,dim(seqs)[2]]), decreasing = TRUE), ], 
	        file = "scoresFor2InputSequences.xls",
		sep = "\t", row.names = FALSE, col.names=TRUE)
	}
	print("Done!")
        if (removegRNADetails[1])
            unlink(outputDir1, recursive = TRUE)
        if (removegRNADetails[2])
            unlink(outputDir2, recursive = TRUE)
	setwd(originalDir)
#scores
	seqs
}
