offTargetAnalysisWithoutBSgenome <-
    function(inputFilePath, format = "fasta", header=FALSE, 
        gRNAoutputName, findgRNAs = TRUE,
        exportAllgRNAs = c("all", "fasta", "genbank", "no"),
        findgRNAsWithREcutOnly = FALSE, 
	REpatternFile = system.file("extdata", "NEBenzymes.fa", 
            package = "CRISPRseek"), 
	minREpatternSize = 4,
	overlap.gRNA.positions = c(17, 18), findPairedgRNAOnly = FALSE, 
        annotatePaired = TRUE, paired.orientation = c("PAMout","PAMin"),
        enable.multicore = FALSE, n.cores.max = 6,
        min.gap = 0, max.gap = 20, gRNA.name.prefix = "",
	PAM.size = 3, gRNA.size = 20, PAM = "NGG", BSgenomeName, 
        chromToSearch = "all", 
        chromToExclude = c("chr17_ctg5_hap1","chr4_ctg9_hap1", "chr6_apd_hap1",
"chr6_cox_hap2", "chr6_dbb_hap3", "chr6_mann_hap4", "chr6_mcf_hap5","chr6_qbl_hap6",
"chr6_ssto_hap7"),
	max.mismatch = 3, 
        PAM.pattern = "NNG$|NGN$", allowed.mismatch.PAM = 1,
        gRNA.pattern = "",
        baseEditing = FALSE, targetBase = "C", editingWindow = 4:8, 
        editingWindow.offtargets = 4:8,
        primeEditing = FALSE, 
        PBS.length = 13L,
	RT.template.length = 8:28,
        RT.template.pattern = "D$",
        corrected.seq,
        targeted.seq.length.change,
        bp.after.target.end = 15L,
        target.start,
        target.end,
        primeEditingPaired.output = "pairedgRNAsForPE.xls",
        min.score = 0, topN = 1000, 
        topN.OfftargetTotalScore = 10, 
        annotateExon = TRUE, txdb, orgAnn, ignore.strand = TRUE, outputDir,
        fetchSequence = TRUE, upstream = 200, downstream = 200,
	upstream.search = 0, downstream.search = 0,
        weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 
        0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583), 
	baseBeforegRNA = 4, baseAfterPAM = 3,
	featureWeightMatrixFile = system.file("extdata", "DoenchNBT2014.csv", 
		package = "CRISPRseek"),
	useScore = TRUE, useEfficacyFromInputSeq = FALSE, 
	outputUniqueREs = TRUE, foldgRNAs = FALSE, 
        gRNA.backbone="GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU",
        temperature = 37,
        overwrite = FALSE,
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
     rule.set = c("Root_RuleSet1_2014", "Root_RuleSet2_2016", "CRISPRscan"),
     mismatch.activity.file = system.file("extdata", 
         "NatureBiot2016SuppTable19DoenchRoot.csv", 
         package = "CRISPRseek"),
     useBSgenome = FALSE, genomeSeqFile,
     predIndelFreq = FALSE,
     predictIndelFreq.onTargetOnly = TRUE,
     method.indelFreq = "Lindel",
     baseBeforegRNA.indelFreq = 13,
     baseAfterPAM.indelFreq = 24
)
{
    cat("Validating input ...\n")
    scoring.method <- match.arg(scoring.method)
    exportAllgRNAs <- match.arg(exportAllgRNAs)
    rule.set <- match.arg(rule.set)
    if (scoring.method ==  "CFDscore") 
    {
        mismatch.activity <- read.csv(mismatch.activity.file)
        required.col <- c("Mismatch.Type", "Position", "Percent.Active")
        if (length(intersect(colnames(mismatch.activity), required.col)) != 
            length(required.col))  
           stop("Please rename the mismatch activity file column to contain at least
              these 3 column names: Mismatch.Type, Position, Percent.Active\n")
    } 
    else if (scoring.method == "Hsu-Zhang")
    {
         if (length(weights) !=  gRNA.size)
             stop("Please make sure the size of weights vector 
                 equals to the gRNA.size!\n")
    }
    if(findgRNAsWithREcutOnly && findgRNAs && !file.exists(REpatternFile))
    {
        stop("Please specify an REpattern file as fasta file with 
            restriction enzyme recognition sequences!")
    }
    if (missing(inputFilePath)) {
        stop("inputFilePath containing the searching sequence, coordinate or a DNAStringSet
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
    if (annotatePaired || findPairedgRNAOnly)
       pairOutputFile <- paste(outputDir, "pairedgRNAs.xls", sep = "")
    REcutDetailFile <- paste(outputDir, "REcutDetails.xls", sep = "")
    bedFile<- paste(outputDir, "gRNAsCRISPRseek.bed", sep = "")
    if (missing(gRNAoutputName) && class(inputFilePath) == "DNAStringSet")
	    stop("Please enter a name for the gRNA ouput file name when DNAStringSet instead of file path provided!")
    if (class(inputFilePath) != "DNAStringSet" && missing(gRNAoutputName))
	    gRNAoutputName = strsplit(basename(inputFilePath), split=".", 
		    fixed=TRUE)[[1]][1]
    if (format =="bed")
    {
        if (missing(BSgenomeName) || class(BSgenomeName) != "BSgenome") {
            stop("BSgenomeName is required as BSgenome object when input file is in bed format!")
        }
        inputFilePath <- getSeqFromBed(inputFilePath, header = header, BSgenomeName = BSgenomeName)
        #### format for filtergRNAs
        format <- "fasta"
    }
    if (findgRNAs)
    {
        cat("Searching for gRNAs ...\n")
	efficacyFile <- paste(outputDir, "gRNAefficacy.xls", sep = "")
	if (chromToSearch == "" || useEfficacyFromInputSeq)
            potential.gRNAs <- findgRNAs(inputFilePath,
               overlap.gRNA.positions = overlap.gRNA.positions,
               baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow,
               primeEditing = primeEditing,
               findPairedgRNAOnly = findPairedgRNAOnly,
               annotatePaired = annotatePaired,
               paired.orientation = paired.orientation,
               pairOutputFile = pairOutputFile, PAM = PAM,
               PAM.location = PAM.location,
               gRNA.pattern = gRNA.pattern, PAM.size = PAM.size,
               gRNA.size = gRNA.size, min.gap = min.gap,
               max.gap = max.gap, name.prefix = gRNA.name.prefix,
               format = format, featureWeightMatrixFile = featureWeightMatrixFile, 
               baseBeforegRNA = baseBeforegRNA, 
	       baseAfterPAM = baseAfterPAM ,
    	       calculategRNAEfficacy = TRUE, efficacyFile = efficacyFile,
               rule.set = rule.set)
         else
	    potential.gRNAs <- findgRNAs(inputFilePath,
               overlap.gRNA.positions = overlap.gRNA.positions,
              baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow,
              primeEditing = primeEditing,
              PBS.length = PBS.length,
                RT.template.length = RT.template.length,
                RT.template.pattern = RT.template.pattern,
                corrected.seq = corrected.seq,
                targeted.seq.length.change = targeted.seq.length.change,
                bp.after.target.end = bp.after.target.end,
                target.start = target.start,
                target.end = target.end,
             primeEditingPaired.output =  primeEditingPaired.output,
               findPairedgRNAOnly = findPairedgRNAOnly,
               annotatePaired = annotatePaired,
               paired.orientation = paired.orientation,
               enable.multicore = enable.multicore,
               n.cores.max = n.cores.max,
               pairOutputFile = pairOutputFile, PAM = PAM,
	       gRNA.pattern = gRNA.pattern, PAM.size = PAM.size,
               PAM.location = PAM.location,
               gRNA.size = gRNA.size, min.gap = min.gap,
               max.gap = max.gap, name.prefix = gRNA.name.prefix, 
               format = format,  rule.set = rule.set)
	if (length(potential.gRNAs) == 0)
        {
		return(cat("no gRNAs found!"))
        }
        
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
                names(subjects) <- gsub( "\t", "", names(subjects))
                names(subjects) <- gsub( "\n", "", names(subjects))
                names(subjects) <- gsub( " ", "", names(subjects))
   	        locuses <- names(subjects)
                
		names.gRNA <- names(potential.gRNAs)
		for (i in 1:length(locuses))
		{
			thisLocus <- gsub("'", "", locuses[i])
        		thisLocus <- gsub(" ", "", thisLocus)
			thisSeq <- tolower(as.character(subjects[[i]]))
			n.bp <- nchar(thisSeq)
			temp <- strsplit(names.gRNA, split=paste(
	         		thisLocus,"_gR",sep=""))
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
					if  (length(grep("f", temp[[j]])) >0)
					{
						temp1 <-strsplit(temp[[j]], "f")
						isForward <- TRUE
					}
					else
					{
						temp1 <-strsplit(temp[[j]], "r")
						isForward <- FALSE
					}
					feature <- temp1[[2]][2]
					feature[is.na(feature)] <- ""
					location <- temp1[[2]][1] 
					if (isForward)
					{
				       	    Start <- location
					    End <- as.numeric(Start) + max(overlap.gRNA.positions) - 
						min(overlap.gRNA.positions)
					    write(paste("     misc_bind       ", Start, "..",
                                                End, sep = ""), append = TRUE, sep="\n",
                                                file = thisFile)
					     write(paste("                     /note=\"gRNAf",
						as.character(feature),
                                                "\"", sep = ""), append = TRUE, sep="\n", file = thisFile)
					}	
					else
					{
                                            End <- location
                                            Start <- as.numeric(End) - max(overlap.gRNA.positions) + 
						min(overlap.gRNA.positions)
					    write(paste("     misc_bind       complement(", 
						    Start, "..", End, ")", sep = ""), 
						    append = TRUE, sep="\n", file = thisFile)
					    write(paste("                     /note=\"gRNAr",
						feature, 
                            			"\"", sep = ""), append = TRUE, sep="\n", file = thisFile)
					}
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
        if ( annotatePaired || findPairedgRNAOnly) 
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
            if (format == "fasta" || format == "fastq")
            {
                potential.gRNAs <- readDNAStringSet(inputFilePath, format,
                      use.names = TRUE)
            }
            else
            {
                stop("format needs to be either fasta,fastq or bed!")
            }
        }
        else
        {
            potential.gRNAs <- inputFilePath
            if (length(names(potential.gRNAs)) == 0)
               names(potential.gRNAs) <- paste("gRNAs", 1:length(potential.gRNAs), sep="")
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
    if (length(chromToSearch) == 1 && chromToSearch == "")
    {
	cat("Done. Please check output files in directory ", outputDir, "\n")
        return(gRNAs)	
    }
    if (useBSgenome && (missing(BSgenomeName) || class(BSgenomeName) != "BSgenome")) {
        stop("BSgenomeName is required as BSgenome object when useBSgenome is set to true for offtarget search!")
    }
    if (!useBSgenome && missing(genomeSeqFile))
        stop("To perform offtarget search with useBSgenome set to FALSE, genomeSeqFile needs to be set to a genome sequence file in fasta format!\n")
    if (annotateExon && (missing(txdb) || (class(txdb) != "TxDb" && 
        class(txdb) != "TranscriptDb")))
    {
        stop("To indicate whether an offtarget is inside an exon, txdb is
            required as TxDb object!")
    }
    names(gRNAs) <- gsub( "\t", "", names(gRNAs))
    names(gRNAs) <- gsub( "\n", "", names(gRNAs))
    names(gRNAs) <- gsub( " ", "", names(gRNAs))

    if (useBSgenome) {
      hits <- searchHits2(gRNAs = gRNAs, PAM = PAM, PAM.pattern = PAM.pattern, 
        BSgenomeName = BSgenomeName, chromToSearch = chromToSearch,
	chromToExclude = chromToExclude,
        max.mismatch = max.mismatch, PAM.size = PAM.size, 
        gRNA.size = gRNA.size, allowed.mismatch.PAM = allowed.mismatch.PAM,
        PAM.location = PAM.location,
        baseEditing = baseEditing, targetBase = targetBase, 
        editingWindow = editingWindow.offtargets) 
    }
    else {
      genomeSeq <- readDNAStringSet(genomeSeqFile)
      if (length(chromToSearch) == 1 && tolower(chromToSearch) == "all")
	   chromInd <- 1:length(genomeSeq)
      else
           chromInd <-  which(names(genomeSeq) %in% chromToSearch)
      chromInd <-  setdiff(chromInd, which(names(genomeSeq) %in% chromToExclude))
      outfile <- tempfile(tmpdir = getwd())
      for (j in chromInd) {
           if (j == chromInd[1]) 
              hits <- searchHits(gRNAs = gRNAs,
                  seqs = genomeSeq[[j]],
                  seqname = names(genomeSeq)[j],
                  max.mismatch = max.mismatch, PAM.size = PAM.size,
                  gRNA.size = gRNA.size, 
                  allowed.mismatch.PAM = allowed.mismatch.PAM,
                  PAM.location = PAM.location,
                  baseEditing = baseEditing, targetBase = targetBase,
                  editingWindow = editingWindow.offtargets,
                  outfile = outfile)
           else
              hits <- rbind(hits, searchHits(gRNAs = gRNAs,
                  PAM = PAM, PAM.pattern = PAM.pattern,
                  seqs = genomeSeq[[j]], 
                  seqname = names(genomeSeq)[j],
                  max.mismatch = max.mismatch, PAM.size = PAM.size,
                  gRNA.size = gRNA.size, 
                  allowed.mismatch.PAM = allowed.mismatch.PAM,
                  PAM.location = PAM.location,
                  baseEditing = baseEditing, targetBase = targetBase,
                  editingWindow = editingWindow.offtargets,
                  outfile = outfile))
      }
   }
if (dim(hits)[1] > 0)
{
    cat("Building feature vectors for scoring ...\n")
    #save(hits, file = "hits.RData")
    featureVectors <- buildFeatureVectorForScoring(hits = hits, 
        canonical.PAM = PAM, gRNA.size = gRNA.size, 
        subPAM.position = subPAM.position, 
        PAM.location = PAM.location, PAM.size = PAM.size)
    cat("Calculating scores ...\n")
    if ( scoring.method ==  "CFDscore")
        scores <- getOfftargetScore2(featureVectors, 
            subPAM.activity = subPAM.activity,
            mismatch.activity.file = mismatch.activity.file)
    else
        scores <- getOfftargetScore(featureVectors, weights = weights)
    #write.table(scores, file="testScore2.xls", sep="\t", row.names=FALSE)
    cat("Annotating, filtering and generating reports ...\n")
    #saveRDS(scores, file="scores.RDS")
    offTargets <- filterOffTargetWithoutBSgenome(scores = scores, outputDir = outputDir,
        BSgenomeName = BSgenomeName, fetchSequence = fetchSequence, txdb = txdb,
            orgAnn = orgAnn, ignore.strand = ignore.strand,
	    min.score = min.score, topN = topN, 
            topN.OfftargetTotalScore = topN.OfftargetTotalScore, 
            upstream = upstream, downstream = downstream, 
            annotateExon = annotateExon, baseBeforegRNA = baseBeforegRNA, 
	    baseAfterPAM = baseAfterPAM, featureWeightMatrixFile = featureWeightMatrixFile,
            rule.set = rule.set, useBSgenome = useBSgenome, genomeSeq = genomeSeq)
#  saveRDS(offTargets, file = "offTargets.RDS") 
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
    if (findgRNAs && (annotatePaired || findPairedgRNAOnly))
    {
        cat("Add paired information...\n")
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
	cat("write gRNAs to bed file...\n")
	on.target <- offTargets$offtargets
	on.target <- unique(subset(on.target, 
            on.target$n.mismatch == 0 & on.target$isCanonicalPAM ==1)) 
	#   as.character(on.target$gRNAPlusPAM) == as.character(on.target$OffTargetSequence)))
	if (dim(on.target)[1] >0)
        {
	   gRNA.bed <- unique(cbind(as.character(on.target$chrom),as.character(on.target$chromStart),
		as.character(on.target$chromEnd), as.character(on.target$name), 
		as.numeric(as.character(on.target$gRNAefficacy)) * 1000, 
		as.character(on.target$strand), 
		as.character(on.target$chromStart), 
		as.character(on.target$chromEnd)))
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
			as.character(on.target$forViewInUCSC),
			as.character(on.target$extendedSequence), 
			as.character(on.target$gRNAefficacy)
                        )) 
	  colnames(on.target) = c("names", "forViewInUCSC", "extendedSequence", "gRNAefficacy")
	  if (useEfficacyFromInputSeq)
	  {
		on.target <- as.data.frame(on.target[,1:2])
		inputEfficacy <- read.table(efficacyFile, sep="\t", header = TRUE, 
			stringsAsFactors=FALSE)
		inputEfficacy <- as.data.frame(cbind(name = inputEfficacy$name,
		        extendedSequence = inputEfficacy$extendedSequence,
			gRNAefficacy = inputEfficacy$gRNAefficacy))
		on.target <- merge(on.target, inputEfficacy, by.x="names", by.y ="name")
	  }
          if(dim(on.target)[1] >0)
	     summary <- unique(merge(on.target, summary, by="names", all=TRUE))
	  write.table(summary[order(as.character(summary$names)), ],
             file = paste(outputDir, "Summary.xls", sep = ""),
             sep = "\t", row.names = FALSE)
	  cat("Scan for REsites in flanking region...\n")
	  if (outputUniqueREs && !missing(BSgenomeName) && 
               class(BSgenomeName) == "BSgenome")
	  {
	    REs.isUnique100 <- uniqueREs(REcutDetails = REcutDetails, 
		   summary = summary, offTargets$offtargets, scanUpstream = 100,
		   scanDownstream =100, BSgenomeName = BSgenomeName)
	    REs.isUnique50 <- uniqueREs(REcutDetails = REcutDetails, 
		   summary = summary, offTargets$offtargets, scanUpstream = 50,
		   scanDownstream = 50, BSgenomeName = BSgenomeName)
	    summary <- cbind(summary, uniqREin200 = REs.isUnique100,
                uniqREin100 = REs.isUnique50) 
            summary$uniqREin200 <- as.character(summary$uniqREin200)
            summary$uniqREin100 <- as.character(summary$uniqREin100)
	  }
	 else
	{
	   REs.isUnique100 = ""
       	   REs.isUnique50 = "" 
	}
    }
    else
    {
       warnings("No on-target found for the input gRNAs with your search criteria!")
       gRNA.bed = ""
       REs.isUnique100 = ""
       REs.isUnique50 = ""
    } 
    if (foldgRNAs)
     {
        source(system.file("extdata/foldgRNAs.R",package = "CRISPRseek"))
	gRNAs.withoutPAM <- substr(as.character(summary$gRNAsPlusPAM), 1, gRNA.size)
        folded.gRNAs <- foldgRNAs(gRNAs.withoutPAM, gRNA.backbone = gRNA.backbone, 
           temperature = temperature)
	if (length(dim(folded.gRNAs)) > 0)
	{
	   if (dim(folded.gRNAs)[1] >1)
	      summary <- cbind(summary, folded.gRNAs[,-1])
	   else
	      summary <- data.frame(c(summary, folded.gRNAs[,-1]))	
	}
     }
    #write.table(summary[order(as.character(summary$forViewInUCSC)), ], 
    ### even there is no perfect target for a gRNA, it will be kept in the summary file
    ### need to calculate the topN offtarget score and distance correctly yet if include those gRNAs without target
     
     gRNAs.notInGenome <- setdiff(names(gRNAs), summary$names)
     if (length(gRNAs.notInGenome) > 0)
     {
         dat2 <- data.frame(matrix(nrow = length(gRNAs.notInGenome), ncol = dim(summary)[2]))
         colnames(dat2) <- colnames(summary)
         dat2$names <- gRNAs.notInGenome
         
         dat2$gRNAsPlusPAM <- paste(substr(as.character(gRNAs[names(gRNAs) %in% gRNAs.notInGenome]), 1, gRNA.size), PAM, sep ="")
         summary <- rbind(summary, dat2)
     }
     if (dim(on.target)[1] == 0)
        write.table(summary[order(as.character(summary$names)), ],
           file = paste(outputDir, "Summary.xls", sep = ""),
           sep = "\t", row.names = FALSE)
     else
        write.table(summary[order(as.character(summary$forViewInUCSC)), ],
           file = paste(outputDir, "Summary.xls", sep = ""), 
           sep = "\t", row.names = FALSE)
    if (predIndelFreq) {
        if (predictIndelFreq.onTargetOnly)
                targets <- unique(subset(offTargets$offtargets,
                     offTargets$offtargets$n.mismatch == 0 & offTargets$offtargets$isCanonicalPAM ==1))
        else
                targets <- offTargets$offtargets
        extendedSequence <- getExtendedSequence(targets,
                 BSgenomeName = BSgenomeName,
                 baseBeforegRNA =  baseBeforegRNA.indelFreq,
                 baseAfterPAM = baseAfterPAM.indelFreq, forMethod = method.indelFreq)
        indelFreq <- predictRelativeFreqIndels(extendedSequence, method = method.indelFreq)
        names(indelFreq) <- paste(targets[,1], targets[,2], targets[,3],
                 targets[,7], targets[,8], targets[,13], round(targets[,19], 3), sep= ",")

        cat("Done. Please check output files in directory \n", outputDir, "\n")
        list(on.target=on.target, summary=summary, offtarget = offTargets$offtargets,
                 gRNAs.bedFormat=gRNA.bed, REcutDetails = REcutDetails,
                 REs.isUnique100 = REs.isUnique100, REs.isUnique50 = REs.isUnique50,
                 indelFreq = indelFreq)
    }
    else { 
        cat("Done. Please check output files in directory \n", outputDir, "\n")
        list(on.target=on.target, summary=summary, offtarget = offTargets$offtargets, 
		 gRNAs.bedFormat=gRNA.bed, REcutDetails = REcutDetails,
		 REs.isUnique100 = REs.isUnique100, REs.isUnique50 = REs.isUnique50)
     }
}
else
{
  x <- paste(substr(as.character(gRNAs), 1, gRNA.size), PAM, sep ="")

  summary <- cbind(names = names(gRNAs), gRNAsPlusPAM = x,top5OfftargetTotalScore = rep("NA", length(gRNAs)),
   	top10OfftargetTotalScore =  rep("NA", length(gRNAs)),
	top1Hit.onTarget.MMdistance2PAM =  rep("NA", length(gRNAs))
      ) 
  write.table(summary,  file = paste(outputDir, "Summary.xls", sep = ""),
        sep = "\t", row.names = FALSE) 
  summary
}
}
