offTargetAnalysis <-
    function(inputFilePath, format = "fasta", findgRNAs = TRUE,
		exportAllgRNAs = c("no", "fasta", "genbank", "all"),
        findgRNAsWithREcutOnly = TRUE, REpatternFile, minREpatternSize = 6,
	overlap.gRNA.positions = c(17, 18), findPairedgRNAOnly = TRUE, 
        min.gap = 0, max.gap = 20, gRNA.name.prefix = "gRNA",
	PAM.size = 3, gRNA.size = 20, PAM = "NGG", BSgenomeName, 
        chromToSearch = "all", max.mismatch = 4, PAM.pattern = "N[A|G]G$",
       gRNA.pattern = "", min.score = 0.5, topN = 100, 
	    topN.OfftargetTotalScore = 10, 
        annotateExon = TRUE, txdb, outputDir,
        fetchSequence = TRUE, upstream = 200, downstream = 200,
        weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 
        0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583), 
        overwrite = FALSE)
{
    cat("Validating input ...\n")
    if(findgRNAsWithREcutOnly && !file.exists(REpatternFile))
    {
        stop("Please specify an REpattern file as fasta file with 
            restriction enzyme recognition sequences!")
    }
    if (missing(inputFilePath)) {
        stop("inputFilePath containing the searching sequence is required!")
    }
    if (!file.exists(inputFilePath)){
        stop("inputfile specified as ", inputFilePath, " does not exists!")
    }
    if (format != "fasta" && format != "fastq") {
        stop("format needs to be either fasta or fastq!")
    }	
    if (substr(outputDir, nchar(outputDir), nchar(outputDir)) != .Platform$file.sep)
    {
        outputDir <- paste(outputDir, "", sep = .Platform$file.sep)
    }
    if (file.exists(outputDir) && ! overwrite)
    {
        cat(outputDir, "exists already. Please type 1 if you want to 
            overwrite the outputDir and 2 if you want to exit.", fill = TRUE)
	input <- readline()
	if(input != 1) { stop("Please change the outputDir!") }
    }
    if (!file.exists(outputDir))
    {
        dir.create(outputDir)
    }
    pairOutputFile <- paste(outputDir, "pairedgRNAs.xls", sep = "")
    REcutDetailFile <- paste(outputDir, "REcutDetails.xls", sep = "")
    if (findgRNAs)
    {
        cat("Searching for gRNAs ...\n")
	potential.gRNAs <- findgRNAs(inputFilePath, 
            findPairedgRNAOnly = findPairedgRNAOnly,
            pairOutputFile = pairOutputFile, PAM = PAM,
	        gRNA.pattern = gRNA.pattern, PAM.size = PAM.size,
            gRNA.size = gRNA.size, min.gap = min.gap, 
            max.gap = max.gap, name.prefix = gRNA.name.prefix)
		if (exportAllgRNAs == "fasta" || exportAllgRNAs == "all")
		{
			writeXStringSet(potential.gRNAs, filepath="allgRNAs.fa")
		}
		if (exportAllgRNAs == "genbank" || exportAllgRNAs == "all")
		{
#### write to genbank format
		}		
	if (findPairedgRNAOnly)
	{
	    gRNAs.RE <- filtergRNAs(potential.gRNAs, 
                pairOutputFile = pairOutputFile, 
                findgRNAsWithREcutOnly = findgRNAsWithREcutOnly,
	        REpatternFile = REpatternFile, 
                format = format,  minREpatternSize = minREpatternSize, 
                overlap.gRNA.positions = overlap.gRNA.positions)
            REcutDetails  <- gRNAs.RE$gRNAREcutDetails
	    write.table(REcutDetails[order(as.character(
                REcutDetails$ForwardgRNAName)), ], file = REcutDetailFile, 
                sep = "\t", row.names = FALSE)		
        }
        else
	{
            gRNAs.RE <- filtergRNAs(potential.gRNAs, 
	        findgRNAsWithREcutOnly = findgRNAsWithREcutOnly,
                REpatternFile = REpatternFile, format = format, 
                minREpatternSize = minREpatternSize, 
                overlap.gRNA.positions = overlap.gRNA.positions)
	    REcutDetails  <- gRNAs.RE$gRNAREcutDetails
	    write.table(REcutDetails[order(as.character(
                REcutDetails$REcutgRNAName)), ], file = REcutDetailFile, 
                sep = "\t", row.names = FALSE)
	}
	if (findgRNAsWithREcutOnly)
	{
	    gRNAs  <- gRNAs.RE$gRNAs
        }
	else
	{
	    gRNAs <- potential.gRNAs
	}
	pairedInformation <- read.table(pairOutputFile, sep = "\t", 
            header = TRUE, stringsAsFactors = FALSE)
    }
    else
    {
        potential.gRNAs <- readDNAStringSet(inputFilePath, format, 
            use.names = TRUE)
	gRNAs.RE <- filtergRNAs(potential.gRNAs, 
            REpatternFile = REpatternFile, format = format, 
            minREpatternSize = minREpatternSize, 
            overlap.gRNA.positions = overlap.gRNA.positions)
	REcutDetails  <- gRNAs.RE$gRNAREcutDetails
	write.table(
            REcutDetails[order(as.character(REcutDetails$REcutgRNAName)), ], 
            file = REcutDetailFile, sep = "\t", row.names = FALSE)
	if (findgRNAsWithREcutOnly)
	{
	    gRNAs  <- gRNAs.RE$gRNAs
	}
	else
	{
	    gRNAs <- potential.gRNAs
	}
	pairedInformation <- ""
    }
    if (chromToSearch == "")
    {
	cat("Done. Please check output files in directory ", outputDir, "\n")
        return(gRNAs)	
    }
    if (missing(BSgenomeName) || class(BSgenomeName) != "BSgenome") {
        stop("BSgenomeName is required as BSgenome object!")
    }
    if (annotateExon && (missing(txdb) || class(txdb) != "TranscriptDb"))
    {
        stop("To indicate whether an offtarget is inside an exon, txdb is
            required as TranscriptDb object!")
    }
    hits <- searchHits(gRNAs = gRNAs, PAM = PAM.pattern, 
        BSgenomeName = BSgenomeName, chromToSearch = chromToSearch, 
        max.mismatch = max.mismatch, PAM.size = PAM.size, 
        gRNA.size = gRNA.size) 
    cat("Building feature vectors for scoring ...\n")
    featureVectors <- buildFeatureVectorForScoring(hits = hits, 
        canonical.PAM = PAM, gRNA.size = gRNA.size)
    cat("Calculating scores ...\n")
    scores <- getOfftargetScore(featureVectors, weights = weights)
    cat("Annotating, filtering and generating reports ...\n")
    offTargets <- filterOffTarget(scores = scores, outputDir = outputDir,
        BSgenomeName = BSgenomeName, fetchSequence = fetchSequence, txdb = txdb,
            min.score = min.score, topN = topN, 
            topN.OfftargetTotalScore = topN.OfftargetTotalScore, 
            upstream = upstream, downstream = downstream, 
            annotateExon = annotateExon)
    summary <- read.table(paste(outputDir, "Summary.xls", sep = ""), sep = "\t", 
        header = TRUE, stringsAsFactors = FALSE) 
    for (i in grep("topOfftarget", names(summary)))
    {
        y <- as.character(summary[,i])
        y[is.na(y)] <- ""
	summary[, i] = y	
    }
    if (findgRNAs)
    {
        PairedgRNAName <- unlist(lapply(1:dim(summary)[1], function(i) {
            as.character(gsub("^\\s+|\\s+$", "", 
                paste(unique(pairedInformation[as.character(
                pairedInformation$ForwardgRNAName) == as.character(
                summary$names[i]),]$ReversegRNAName),
                unique(pairedInformation[as.character(
                pairedInformation$ReversegRNAName) == as.character(
                summary$names[i]),]$ForwardgRNAName),
                collapse = " ")))
        }))
    }
    if (findPairedgRNAOnly && findgRNAs)
    {
        REname <- unlist(lapply(1:dim(summary)[1], function(i) {
            gsub("^\\s+|\\s+$", "", gsub("NA", "", 
                paste(unique(REcutDetails[as.character(
                REcutDetails$ForwardREcutgRNAName) == as.character(
                summary$names[i]),]$ForwardREname),
                unique(REcutDetails[as.character(
                REcutDetails$ReverseREcutgRNAName) == 
                as.character(summary$names[i]), ]$ReverseREname), 
                collapse = " ")))
       }))
       summary <- cbind(summary, PairedgRNAName, REname)
    }
    else
    {
        REname <- unlist(lapply(1:dim(summary)[1], function(i) {
            gsub("^\\s+|\\s+$", "", gsub("NA", "", paste(unique(
                REcutDetails[as.character(REcutDetails$REcutgRNAName) == 
                as.character(summary$names[i]), ]$REname), collapse = " ")))
        }))
        summary <- cbind(summary, REname)
    }
    write.table(summary[order(as.character(summary$names)), ], 
        file = paste(outputDir, "Summary.xls", sep = ""), 
        sep = "\t", row.names = FALSE)
    cat("Done. Please check output files in directory ", outputDir, "\n")
}
