offTargetAnalysis <-
    function(inputFilePath, format = "fasta", gRNAoutputName, findgRNAs = TRUE,
		exportAllgRNAs = c("all", "fasta", "genbank", "no"),
        findgRNAsWithREcutOnly = TRUE, 
	REpatternFile = system.file("extdata", "NEBenzymes.fa", 
            package = "CRISPRseek"), 
	minREpatternSize = 6,
	overlap.gRNA.positions = c(17, 18), findPairedgRNAOnly = TRUE, 
        min.gap = 0, max.gap = 20, gRNA.name.prefix = "gRNA",
	PAM.size = 3, gRNA.size = 20, PAM = "NGG", BSgenomeName, 
        chromToSearch = "all", max.mismatch = 3, PAM.pattern = "N[A|G]G$",
        gRNA.pattern = "", min.score = 0.5, topN = 100, 
        topN.OfftargetTotalScore = 10, 
        annotateExon = TRUE, txdb, orgAnn, outputDir,
        fetchSequence = TRUE, upstream = 200, downstream = 200,
        weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 
        0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583), 
	baseBeforegRNA = 4, baseAfterPAM = 3,
	featureWeightMatrixFile = system.file("extdata", "DoenchNBT2014.csv", 
		package = "CRISPRseek"),
	useScore = TRUE,
        overwrite = FALSE)
{
    cat("Validating input ...\n")
    if(findgRNAsWithREcutOnly && findgRNAs && !file.exists(REpatternFile))
    {
        stop("Please specify an REpattern file as fasta file with 
            restriction enzyme recognition sequences!")
    }
    if (missing(inputFilePath)) {
        stop("inputFilePath containing the searching sequence or a DNAStringSet
             object is required!")
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
    bedFile<- paste(outputDir, "gRNAsCRISPRseek.bed", sep = "")
    if (missing(gRNAoutputName) && class(inputFilePath) == "DNAStringSet")
	stop("Please enter a name for the gRNA ouput file name!")
    if (class(inputFilePath) != "DNAStringSet" && missing(gRNAoutputName))
	gRNAoutputName = strsplit(basename(inputFilePath), split=".", 
		fixed=TRUE)[[1]][1]
    if (findgRNAs)
    {
        cat("Searching for gRNAs ...\n")
	potential.gRNAs <- findgRNAs(inputFilePath, 
            findPairedgRNAOnly = findPairedgRNAOnly,
            pairOutputFile = pairOutputFile, PAM = PAM,
	    gRNA.pattern = gRNA.pattern, PAM.size = PAM.size,
            gRNA.size = gRNA.size, min.gap = min.gap, 
            max.gap = max.gap, name.prefix = gRNA.name.prefix, format = format)
	if (length(potential.gRNAs) == 0)
		stop("no gRNAs found!")
	if (length(potential.gRNAs) > 0 && (exportAllgRNAs == "fasta" || exportAllgRNAs == "all"))
	{
		writeXStringSet(potential.gRNAs, filepath= file.path(outputDir,
                     paste(gRNAoutputName,"allgRNAs.fa", sep="")))
	}
	if (length(potential.gRNAs) > 0 && (exportAllgRNAs == "genbank" || exportAllgRNAs == "all"))
	{
		if (class(inputFilePath) == "DNAStringSet")
			subjects <- inputFilePath
		else
			subjects <- readDNAStringSet(inputFilePath, format=format,
				use.names = TRUE)
   	        locuses <- names(subjects)
		names.gRNA <- names(potential.gRNAs)
		for (i in 1:length(locuses))
		{
			thisLocus <- gsub("'", "", locuses[i])
        		thisLocus <- gsub(" ", "", thisLocus)
			thisSeq <- tolower(as.character(subjects[[i]]))
			n.bp <- nchar(thisSeq)
			temp <- strsplit(names.gRNA, split=paste("_", 
	         		thisLocus,"Start",sep=""))
			locus <- paste("LOCUS       ", thisLocus,
                                        "                     ", n.bp,
					" bp    dna     linear   UNK", sep="")
			definition <- paste("DEFINITION  CRISPRseek output for ",
				    gRNAoutputName, " sequence", sep = "")
			accession <- "ACCESSION   unknown"
			features <- "FEATURES             Location/Qualifiers"
			header = rbind(locus, definition, accession, features)
			found.gRNA <- 0
			for (j in 1:length(temp))
			{
				if (length(temp[[j]]) >1){
					found.gRNA <- found.gRNA + 1
					if (found.gRNA == 1)
					{
					    thisFile <- file.path(outputDir,
					      paste(thisLocus, "gbk", sep="."))
                            		    write(header, thisFile)
					}
					feature <- temp[[j]][1]
					location <- strsplit(temp[[j]][2], split="End")
                        		Start <- min(as.numeric(location[[1]]))
					End <- max(as.numeric(location[[1]]))
					if (Start == as.numeric(location[[1]])[2])
					    write(paste("     misc_bind       complement(", 
						    (Start + PAM.size), "..", End, ")", sep = ""), 
						    append = TRUE, sep="\n", file = thisFile)
					else	
					    write(paste("     misc_bind       ", Start, "..",
						(End - PAM.size), sep = ""), append = TRUE, sep="\n",
                                		file = thisFile)
					write(paste("                     /note=\"", feature, 
                            			"\"", sep = ""), append = TRUE, sep="\n", file = thisFile)
					}
			}
			if (found.gRNA > 0){
			    write("ORIGIN", append = TRUE, sep="\n", file = thisFile)
                    	    seq.lines <- floor(nchar(thisSeq) / 60) + 1
                            for (k in 1:seq.lines) {
                                line.start <- (k - 1) * 60 + 1
                        	line.end <- min(line.start + 59, nchar(thisSeq))
                        	n.leading.spaces <- 9 - nchar(line.start)
                        	leading.spaces <- paste(rep(" ", n.leading.spaces), 
                            		collapse = "")
                        	seq.thisLine <- substr(thisSeq, line.start, line.end)
                        	len.thisLine <- nchar(seq.thisLine)
                        	n.seg <- floor(len.thisLine /10) + 1
                        	for (l in 1:n.seg) {
                            		seg.start <- (l -1) * 10 + 1
                            		seg.end <- min(seg.start + 9, len.thisLine)
                            		if (l == 1)
                                		seq.thisLine.formatted <- substr(seq.thisLine,
                                    			seg.start, seg.end)
                            		else
                                	seq.thisLine.formatted <- paste(
                                    		seq.thisLine.formatted,
                                    		substr(seq.thisLine, seg.start, seg.end),
                                    		sep = " ")
                             	}
                        	write(paste(leading.spaces, line.start, " ", 
                            		seq.thisLine.formatted, sep = ""),
                            		append = TRUE, sep="\n", file = thisFile)
                    	}
				write("//", append = TRUE, sep="\n", file = thisFile)
			}
		    }
		}
	if (findPairedgRNAOnly && length(potential.gRNAs) >0)
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
        else if (length(potential.gRNAs) >0)
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
        if (class(inputFilePath) != "DNAStringSet")
        {
            if (! file.exists(inputFilePath)) {
                stop("inputfile specified as ", inputFilePath, " does not exists!")
            }
            if (format != "fasta" && format != "fastq") 
            {
                stop("format needs to be either fasta or fastq!")
            }
            potential.gRNAs <- readDNAStringSet(inputFilePath, format, 
                use.names = TRUE)
        }
        else
        {
            potential.gRNAs <- inputFilePath
        }
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
    if (annotateExon && (missing(txdb) || (class(txdb) != "TxDb" && class(txdb) != "TranscriptDb")))
    {
        stop("To indicate whether an offtarget is inside an exon, txdb is
            required as TxDb object!")
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
            orgAnn = orgAnn, min.score = min.score, topN = topN, 
            topN.OfftargetTotalScore = topN.OfftargetTotalScore, 
            upstream = upstream, downstream = downstream, 
            annotateExon = annotateExon, baseBeforegRNA = baseBeforegRNA, 
	    baseAfterPAM = baseAfterPAM, featureWeightMatrixFile = featureWeightMatrixFile)
    cat("Done annotating\n")
    summary <- read.table(paste(outputDir, "Summary.xls", sep = ""), sep = "\t", 
        header = TRUE, stringsAsFactors = FALSE) 
    if (dim(summary)[2] == 1)
    	summary <- as.data.frame(t(data.matrix(offTargets$summary))) 
    for (i in grep("topOfftarget", names(summary)))
    {
        y <- as.character(summary[,i])
        y[is.na(y)] <- ""
	summary[, i] = y	
    }
    cat("Add paired information...\n")
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
    cat("Add RE information...\n")
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
	seq <- as.character(summary$gRNAsPlusPAM)
	#n.C <- unlist(lapply(1:length(seq), function(i) {
	#					 table(factor(s2c(substr(seq[i],1,gRNA.size)), levels=c("C")))
	#					 }))
	#n.G <- unlist(lapply(1:length(seq), function(i) {
	#					 table(factor(s2c(substr(seq[i],1,gRNA.size)), levels=c("G")))
	#					 }))
	#summary <- cbind(summary,GC.percent.gRNA = (n.G+n.C)/gRNA.size * 100)
	#numberOfT.last4Position.gRNA <- unlist(lapply(1:length(seq), function(i) {
	#	 table(factor(s2c(substr(seq[i],gRNA.size-3,gRNA.size)), levels=c("T")))
	#	 }))
	#summary <- cbind(summary, numberOfT.last4Position.gRNA )
	cat("write to gRNAs to bed file...\n")
	on.target <- offTargets$offtargets
	on.target <- unique(subset(on.target, on.target$n.mismatch == 0 & on.target$NGG == 1))
	gRNA.bed <- unique(cbind(as.character(on.target$chrom),on.target$chromStart,
		on.target$chromEnd, as.character(on.target$name), 
		on.target$gRNAefficiency * 1000, 
		on.target$strand, on.target$chromStart, on.target$chromEnd))
	if (!useScore)
	{
		gRNA.bed <- cbind(gRNA.bed, rep("255,0,0",dim(gRNA.bed)[1]))	
		gRNA.bed[gRNA.bed[,6] == "-",9] = "0,255,0"
	}
	#### UCSC genome browser is 0-based instead of 1 based index
	gRNA.bed[, 2] = as.numeric(gRNA.bed[, 2]) -1
	gRNA.bed[, 3] = as.numeric(gRNA.bed[, 3])
	gRNA.bed[gRNA.bed[,6] == "+" ,7] <- as.numeric(gRNA.bed[gRNA.bed[,6] == "+" ,2]) + 
		min(overlap.gRNA.positions) - 1
        gRNA.bed[gRNA.bed[,6] == "-" ,7] <- as.numeric(gRNA.bed[gRNA.bed[,6] == "-" ,3]) - 
		max(overlap.gRNA.positions) 
        gRNA.bed[gRNA.bed[,6] == "+", 8] <- as.numeric(gRNA.bed[gRNA.bed[,6] == "+", 2]) + 
		max(overlap.gRNA.positions)
	gRNA.bed[gRNA.bed[,6] == "-", 8] <- as.numeric(gRNA.bed[gRNA.bed[,6] == "-", 3]) -
		 min(overlap.gRNA.positions) + 1 
	write.table("track name=\"gRNA sites\" description=\"CRISPRseek\" visibility=2 useScore=1 itemRgb=\"On\"", file=bedFile, col.names=FALSE, row.names=FALSE, quote = FALSE)
	write.table(gRNA.bed, file=bedFile, sep=" ", row.names=FALSE, col.names=FALSE, append=TRUE, quote = FALSE)
	on.target <- unique(cbind(as.character(on.target$name),
                           as.character(on.target$toViewInUCSC),
                          as.character(on.target$extendedSequence), on.target$gRNAefficiency
                        )) 
	colnames(on.target) = c("names", "toViewInUCSC", "extendedSequence", "gRNAefficiency")
	summary <- unique(merge(on.target, summary, by="names", all.y=TRUE))
	
    write.table(summary[order(as.character(summary$names)), ], 
        file = paste(outputDir, "Summary.xls", sep = ""), 
        sep = "\t", row.names = FALSE)
    cat("Done. Please check output files in directory ", outputDir, "\n")
    list(summary=summary, offtarget=offTargets$offtargets, gRNAs=gRNA.bed)
}
