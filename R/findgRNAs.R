findgRNAs <-
    function (inputFilePath, format = "fasta", PAM = "NGG", PAM.size = 3, 
        findPairedgRNAOnly = FALSE, gRNA.pattern = "", gRNA.size = 20, 
        min.gap = 0, max.gap = 20, pairOutputFile, name.prefix = "gRNA",
	featureWeightMatrixFile = system.file("extdata", 
        "DoenchNBT2014.csv", package = "CRISPRseek"), baseBeforegRNA = 4, 
	baseAfterPAM = 3,
	calculategRNAEfficacy = FALSE, efficacyFile) 
{
    if (missing(inputFilePath)) {
        stop("inputFilePath containing the searching sequence 
	   or a DNAStringSet object is required!")
    }
    if (class(inputFilePath) != "DNAStringSet")
    {
    	if (! file.exists(inputFilePath)) {
           stop("inputfile specified as ", inputFilePath, " does not exists!")
	}
    	if (format != "fasta" && format != "fastq") 
    	{
       	   stop("format needs to be either fasta or fastq!")
    	}
    	subjects <- readDNAStringSet(inputFilePath, format, use.names = TRUE)
    }
    else
    {
	subjects <- inputFilePath
     }
    if (length(subjects) == 0)
    {
        stop("The input file contains no sequence! This could be caused by 
            wrong format of the file. If file is created in mac, you could 
            reformat to text by typing tr \"\\r\" \"\\n\" >newfile in the 
            command line")
    }
    #for (i in 1:length(subjects))
    toAppend = FALSE
    colNames = TRUE
    all.gRNAs.df <- do.call(rbind, lapply(1:length(subjects), function(p)
   {
        subjectname <- gsub("'", "", names(subjects)[p])
        subjectname <- gsub(" ", "", subjectname)
        subjectname <- gsub("\t", ":", subjectname)
        subject <- subjects[[p]]
        revsubject <- reverseComplement(subject)
        PAM <- translatePattern(PAM)
        PAM <- paste(PAM, "$", sep = "")
	gRNA.pattern <- translatePattern(gRNA.pattern)
        seq.len <- nchar(as.character(subject))
        x <- seq.len - gRNA.size - PAM.size + 1
        pos.gRNAs <- do.call(rbind, lapply(1:x, function(i){
            seq <- subseq(as.character(subject), i, 
                (i + gRNA.size + PAM.size -1))
	    extendedSequence <- seq
            if (calculategRNAEfficacy)
            {
	        extended.start <- max(1, (i - baseBeforegRNA))
	        extended.end <- i + gRNA.size + PAM.size -1 + baseAfterPAM
                extended.end <- min(extended.end, length(subject))	
	        extendedSequence <- subseq(as.character(subject), extended.start,
                   extended.end)
            }
            if (regexpr(PAM, seq, perl = TRUE,ignore.case = TRUE)[[1]] > 0 &&
		(gRNA.pattern =="" || regexpr(gRNA.pattern, seq, perl = TRUE,
		ignore.case = TRUE)[[1]] > 0))
                c(seq, paste(subjectname,"Start",i,"End", 
                    (i + gRNA.size + PAM.size-1), sep = ""),i,"+", extendedSequence)
        }))
        minus.gRNAs <- do.call(rbind, lapply(1:x, function(i){
            seq <- subseq(as.character(revsubject), i, 
                (i + gRNA.size + PAM.size -1))
	    extendedSequence <- seq
            if (calculategRNAEfficacy)
	    {
	    	extended.start <- max(1, (i - baseBeforegRNA))
            	extended.end <- i + gRNA.size + PAM.size -1 + baseAfterPAM
           	extended.end <- min(extended.end, length(subject))
            	extendedSequence <- subseq(as.character(revsubject), extended.start,
                    extended.end)
	    }
            if (regexpr(PAM, seq, perl = TRUE,ignore.case = TRUE)[[1]] > 0 &&
		(gRNA.pattern =="" || regexpr(gRNA.pattern, seq, perl = TRUE,
	     	    ignore.case = TRUE)[[1]] > 0))
                c(seq, paste( subjectname,"Start", (seq.len - i + 1), "End", 
                    (seq.len - i + 1 - gRNA.size - PAM.size + 1), sep = ""),
                    (seq.len - i + 1), "-", extendedSequence)			
        }))
        plus.index <- numeric()
        minus.index <- numeric()
        forward.index <- 0
        reverse.index <- 0
        if (length(pos.gRNAs) > 1)
            n.plus.gRNAs <- dim(pos.gRNAs)[1]
        else
            n.plus.gRNAs <- 0
        if (length(minus.gRNAs) > 1)
            n.minus.gRNAs <- dim(minus.gRNAs)[1]
        else
            n.minus.gRNAs <- 0
        if (findPairedgRNAOnly && (n.minus.gRNAs  * n.plus.gRNAs) == 0)
        {
            warning(paste("No paired gRNAs found in the input sequence", subjectname))
            all.gRNAs <- DNAStringSet()
        }
	else
	{
            temp = matrix(nrow =0, ncol=5)
            colnames(temp)[1:5] <- c( "ReversegRNAPlusPAM",
                "ReversegRNAName", "ForwardgRNAPlusPAM",
                "ForwardgRNAName", "gap")
	        if (p >1)
	        {
	            toAppend = TRUE
	            colNames = FALSE
	        }
            write.table(temp, file = pairOutputFile, sep = "\t",
                row.names = FALSE, quote = FALSE, append = toAppend, 
	            col.names = colNames)
            if (n.minus.gRNAs > 0 && n.plus.gRNAs > 0)
            {
                for (j in 1:n.plus.gRNAs)
                {
                    for (k in 1:n.minus.gRNAs)
                    {
                        if ((as.numeric(as.character(pos.gRNAs[j,3])) - 1 - 
                            as.numeric(as.character(minus.gRNAs[k,3]))) >= min.gap
                            && (as.numeric(as.character(pos.gRNAs[j,3])) - 
                            as.numeric(as.character(minus.gRNAs[k,3]))) <= 
                            max.gap)
                        {
                            if(findPairedgRNAOnly)
                            {
                                if (!k %in% minus.index)
                                {
                                    reverse.index <- reverse.index + 1
                                    minus.gRNAs[k,2] <- paste(name.prefix, "r", 
                                        reverse.index, "_", minus.gRNAs[k, 2], 
                                        sep = "")
                                }
                                if (!j %in% plus.index)
                                {
                                    forward.index <- forward.index + 1
                                    pos.gRNAs[j,2] <- paste(name.prefix,"f",
                                        forward.index, "_", pos.gRNAs[j, 2], 
                                        sep = "")
                                }
                            }
                            plus.index <- c(plus.index, j)
                            minus.index <- c(minus.index, k)
                        } ### if paired
                    } #for k
                }## for j
                if (length(minus.index) > 0 && length(plus.index) > 0)
                {
                    paired <- cbind(minus.gRNAs[minus.index,1], 
                        minus.gRNAs[minus.index, 2], pos.gRNAs[plus.index,1], 
                        pos.gRNAs[plus.index, 2],
                        gap = as.numeric(as.character(
                        pos.gRNAs[plus.index,3])) - 1 - as.numeric(
                        as.character(minus.gRNAs[minus.index, 3])))
                    colnames(paired)[1:5] <- c( "ReversegRNAPlusPAM", 
                        "ReversegRNAName", "ForwardgRNAPlusPAM", 
                        "ForwardgRNAName", "gap")
                    if(findPairedgRNAOnly)
                    {
                        plus.index <- unique(plus.index)
                        minus.index <- unique(minus.index)
                        all.gRNAs <- DNAStringSet(c(pos.gRNAs[plus.index, 1], 
                            minus.gRNAs[minus.index,1]))
                        names(all.gRNAs) <- c(pos.gRNAs[plus.index,2], 
                            minus.gRNAs[minus.index,2])
                    }
                    else
                    {
                        paired[, 2] <- paste(name.prefix, "r", minus.index, "_", 
                            paired[, 2], sep = "")
                        paired[, 4] <- paste(name.prefix, "f", plus.index, "_", 
                            paired[, 4], sep = "")
                    }
                    if (dim(paired)[1] == 1)
                        write.table(paired, file = pairOutputFile, sep = "\t", 
                            row.names = FALSE, quote = FALSE, append = toAppend,
			                col.names = colNames)
                    else
                        write.table(paired[order(as.character(paired[,4])), ], 
                            file = pairOutputFile, sep = "\t", row.names = FALSE,
                            quote = FALSE, append = toAppend, col.names = colNames)
                }
                else if (findPairedgRNAOnly)
                {
                    warning(paste("No paired gRNAs found for sequence",
                        subjectname))
                   all.gRNAs <- DNAStringSet()
                }
            }### if pos.gRNAs and minus.gRNAs not empty
            if (! findPairedgRNAOnly)
            {
                all.gRNAs <- DNAStringSet(c(pos.gRNAs[,1], minus.gRNAs[,1]))
                if (length(all.gRNAs) == 0)
                    warning(paste("No gRNAs found in the input sequence", 
                        subjectname))	
                if (n.plus.gRNAs == 0 && n.minus.gRNAs >0 )
                    names(all.gRNAs) <- paste(name.prefix,
                        c(rep("r", dim(minus.gRNAs)[1])),
                        c(1:dim(minus.gRNAs)[1]), "_", 
                        c(minus.gRNAs[,2]), sep = "")
                if (n.minus.gRNAs == 0 && n.plus.gRNAs > 0)
                    names(all.gRNAs) <- paste(name.prefix,c(rep("f", 
                        dim(pos.gRNAs)[1])),
                        c(1:dim(pos.gRNAs)[1]), 
                        c(pos.gRNAs[,2]), sep = "")
                if (n.minus.gRNAs > 0 && n.plus.gRNAs > 0)
                    names(all.gRNAs) <- paste(name.prefix,c(rep("f", 
                        dim(pos.gRNAs)[1]), rep("r", dim(minus.gRNAs)[1])),
                        c(1:dim(pos.gRNAs)[1], 1:dim(minus.gRNAs)[1]), "_", 
                        c(pos.gRNAs[,2], minus.gRNAs[,2]), sep = "")
            }
	 } ### no paired found and findPairedOnly   
	 forEffi <- rbind(pos.gRNAs, minus.gRNAs)
	 forEffi <- subset(forEffi, forEffi[,1] %in% as.character(all.gRNAs))	
	 if (length(all.gRNAs) >0)
	     #cbind(as.data.frame(all.gRNAs), names(all.gRNAs), forEffi)
	     forEffi
    })) ## do call subjects
    if (calculategRNAEfficacy)
    {
        featureWeightMatrix <- read.csv(featureWeightMatrixFile, header=TRUE)
        effi <- calculategRNAEfficiency(all.gRNAs.df[,5], 
	    baseBeforegRNA = baseBeforegRNA, 
	    featureWeightMatrix = featureWeightMatrix, gRNA.size = gRNA.size) 
        extendedSequences <- cbind(all.gRNAs.df, effi)
        colnames(extendedSequences)  = c("gRNAplusPAM", "name", "start", "strand", 
	    "extendedSequence", "gRNAefficacy")
        extendedSequences[nchar(extendedSequences[,5]) 
	    < baseBeforegRNA + gRNA.size + PAM.size + baseAfterPAM, 6] <- 
            "extended sequence too short"
	write.table(extendedSequences[order(extendedSequences[,6], decreasing = TRUE),],
             file = efficacyFile, sep="\t", row.names = FALSE)
    }
    #else
    #{
#	extendedSequences <- all.gRNAs.df[,1:4]
#	colnames(extendedSequences)  = c("gRNAplusPAM", "name", "start", "strand")
#    }
    all.gRNAs <- DNAStringSet(all.gRNAs.df[,1])
    names(all.gRNAs) = all.gRNAs.df[,2]
    #list(all.gRNAs=all.gRNAs, extendedSequences = extendedSequences)
    all.gRNAs
}
