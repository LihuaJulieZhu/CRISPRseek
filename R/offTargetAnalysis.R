offTargetAnalysis <-
    function(inputFilePath, format = "fasta", findSpacers = TRUE, 
        findSpacersWithREcutOnly = TRUE, REpatternFile, minREpatternSize = 6,
	overlap.spacer.positions = c(17, 18), findPairedSpacerOnly = TRUE, 
        min.gap = 0, max.gap = 20, spacer.name.prefix = "gRNA",
	PAM.size = 3, spacer.size = 20, PAM = "NGG", BSgenomeName, 
        chromToSearch = "all", max.mismatch = 4, PAM.pattern = "N[A|G]G$",
        min.score = 0.5, topN = 100, topN.OfftargetTotalScore = 10, 
        annotateExon = TRUE, txdb, outputDir,
        fetchSequence = TRUE, upstream = 200, downstream = 200,
        weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 
        0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583), 
        overwrite = FALSE)
{
    cat("Validating input ...\n")
    if(findSpacersWithREcutOnly && !file.exists(REpatternFile))
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
    pairOutputFile <- paste(outputDir, "pairedSpacers.xls", sep = "")
    REcutDetailFile <- paste(outputDir, "REcutDetails.xls", sep = "")
    if (findSpacers)
    {
        cat("Searching for spacers ...\n")
	potential.spacers <- findSpacers(inputFilePath, 
            findPairedSpacerOnly = findPairedSpacerOnly,
            pairOutputFile = pairOutputFile, PAM = PAM, 
            PAM.size = PAM.size, spacer.size = spacer.size, min.gap = min.gap, 
            max.gap = max.gap, name.prefix = spacer.name.prefix)
	if (findPairedSpacerOnly)
	{
	    spacers.RE <- filterSpacers(potential.spacers, 
                pairOutputFile = pairOutputFile, 
                findSpacersWithREcutOnly = findSpacersWithREcutOnly,
	        REpatternFile = REpatternFile, 
                format = format,  minREpatternSize = minREpatternSize, 
                overlap.spacer.positions = overlap.spacer.positions)
            REcutDetails  <- spacers.RE$spacerREcutDetails
	    write.table(REcutDetails[order(as.character(
                REcutDetails$ForwardSpacerName)), ], file = REcutDetailFile, 
                sep = "\t", row.names = FALSE)		
        }
        else
	{
            spacers.RE <- filterSpacers(potential.spacers, 
	        findSpacersWithREcutOnly = findSpacersWithREcutOnly,
                REpatternFile = REpatternFile, format = format, 
                minREpatternSize = minREpatternSize, 
                overlap.spacer.positions = overlap.spacer.positions)
	    REcutDetails  <- spacers.RE$spacerREcutDetails
	    write.table(REcutDetails[order(as.character(
                REcutDetails$REcutSpacerName)), ], file = REcutDetailFile, 
                sep = "\t", row.names = FALSE)
	}
	if (findSpacersWithREcutOnly)
	{
	    spacers  <- spacers.RE$spacers
        }
	else
	{
	    spacers <- potential.spacers
	}
	pairedInformation <- read.table(pairOutputFile, sep = "\t", 
            header = TRUE, stringsAsFactors = FALSE)
    }
    else
    {
        potential.spacers <- readDNAStringSet(inputFilePath, format, 
            use.names = TRUE)
	spacers.RE <- filterSpacers(potential.spacers, 
            REpatternFile = REpatternFile, format = format, 
            minREpatternSize = minREpatternSize, 
            overlap.spacer.positions = overlap.spacer.positions)
	REcutDetails  <- spacers.RE$spacerREcutDetails
	write.table(
            REcutDetails[order(as.character(REcutDetails$REcutSpacerName)), ], 
            file = REcutDetailFile, sep = "\t", row.names = FALSE)
	if (findSpacersWithREcutOnly)
	{
	    spacers  <- spacers.RE$spacers
	}
	else
	{
	    spacers <- potential.spacers
	}
	pairedInformation <- ""
    }
    if (chromToSearch == "")
    {
	cat("Done. Please check output files in directory ", outputDir, "\n")
        return(spacers)	
    }
    if (missing(BSgenomeName) || class(BSgenomeName) != "BSgenome") {
        stop("BSgenomeName is required as BSgenome object!")
    }
    if (annotateExon && (missing(txdb) || class(txdb) != "TranscriptDb"))
    {
        stop("To indicate whether an offtarget is inside an exon, txdb is
            required as TranscriptDb object!")
    }
    hits <- searchHits(spacers = spacers, PAM = PAM.pattern, 
        BSgenomeName = BSgenomeName, chromToSearch = chromToSearch, 
        max.mismatch = max.mismatch, PAM.size = PAM.size, 
        spacer.size = spacer.size) 
    cat("Building feature vectors for scoring ...\n")
    featureVectors <- buildFeatureVectorForScoring(hits = hits, 
        canonical.PAM = PAM, spacer.size = spacer.size)
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
    if (findSpacers)
    {
        PairedSpacerName <- unlist(lapply(1:dim(summary)[1], function(i) {
            as.character(gsub("^\\s+|\\s+$", "", 
                paste(unique(pairedInformation[as.character(
                pairedInformation$ForwardSpacerName) == as.character(
                summary$names[i]),]$ReverseSpacerName),
                unique(pairedInformation[as.character(
                pairedInformation$ReverseSpacerName) == as.character(
                summary$names[i]),]$ForwardSpacerName),
                collapse = " ")))
        }))
    }
    if (findPairedSpacerOnly && findSpacers)
    {
        REname <- unlist(lapply(1:dim(summary)[1], function(i) {
            gsub("^\\s+|\\s+$", "", gsub("NA", "", 
                paste(unique(REcutDetails[as.character(
                REcutDetails$ForwardREcutSpacerName) == as.character(
                summary$names[i]),]$ForwardREname),
                unique(REcutDetails[as.character(
                REcutDetails$ReverseREcutSpacerName) == 
                as.character(summary$names[i]), ]$ReverseREname), 
                collapse = " ")))
       }))
       summary <- cbind(summary, PairedSpacerName, REname)
    }
    else
    {
        REname <- unlist(lapply(1:dim(summary)[1], function(i) {
            gsub("^\\s+|\\s+$", "", gsub("NA", "", paste(unique(
                REcutDetails[as.character(REcutDetails$REcutSpacerName) == 
                as.character(summary$names[i]), ]$REname), collapse = " ")))
        }))
        summary <- cbind(summary, REname)
    }
    write.table(summary[order(as.character(summary$names)), ], 
        file = paste(outputDir, "Summary.xls", sep = ""), 
        sep = "\t", row.names = FALSE)
    cat("Done. Please check output files in directory ", outputDir, "\n")
}
