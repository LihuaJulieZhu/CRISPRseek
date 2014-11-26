filtergRNAs <-
    function(all.gRNAs, pairOutputFile = "", 
        findgRNAsWithREcutOnly = FALSE, 
	REpatternFile = system.file("extdata", "NEBenzymes.fa",
            package = "CRISPRseek"), format = "fasta",
        minREpatternSize = 4, overlap.gRNA.positions = c(17, 18),
		overlap.allpos = TRUE)
{
	if (length(all.gRNAs) == 0)
		stop("all.gRNAs contains no gRNAs!")
    if (missing(REpatternFile)) {
        stop("REpatternFile containing the restriction enzyme cut pattern 
            is required!")
    }
    if (!file.exists(REpatternFile)) {
        stop("REpatternFile specified as ", REpatternFile, " does not exists!")
    }
    if (format != "fasta" && format != "fastq") {
        stop("format needs to be either fasta or fastq!")
    }
    patterns <- unique(readDNAStringSet(REpatternFile, format, 
        use.names = TRUE))
    patterns <- patterns[width(patterns) >= minREpatternSize, ]
    countRE <- as.data.frame(table(names(patterns)))
    singleRE <- patterns[names(patterns) %in% countRE[countRE[,2] == 1, 1]]
    duplicateRE <- patterns[names(patterns) %in% countRE[countRE[,2] > 1, 1]]
    duplicateRE <- duplicateRE[order(names(duplicateRE))]
    if (length(duplicateRE) > 0)
    {
        print("More than one pattern found for the following REs and these RE
            will be skipped! If you want to include them in the search, please 
            correct the REpattern file and run this again!")
        print(unique(names(duplicateRE)))
    }
    patterns <- singleRE
    seqs <- as.character(all.gRNAs)
    seq.names <- names(all.gRNAs)
    min.pStart.plus <- min(overlap.gRNA.positions)
    max.pStart.plus <- max(overlap.gRNA.positions)
    gRNAs.RE <- data.frame(gRNAPlusPAM = "", REcutgRNAName = "", 
        REname = "", REpattern = "", REcutStart = "", REcutEnd = "")
    for (j in 1:length(patterns))
    {
        pattern.name <- gsub("'", "", names(patterns)[j])
        pattern <- patterns[[j]]
		this.pattern.size <- length(pattern)
        revpattern <- reverseComplement(pattern)
        if (revpattern != pattern)
        {
			revpattern <- translatePattern(revpattern)
            minus.gRNAs <- do.call(rbind, lapply(1:length(all.gRNAs), 
                function(i){
                    res1  <- as.numeric(gregexpr(revpattern, seqs[i],
                        perl = TRUE,fixed = FALSE,ignore.case = TRUE)[[1]])
                    do.call(rbind, lapply(1:length(res1), function(k) {
						if (res1[k] >0 && !overlap.allpos && 
							((min.pStart.plus >= res1[k] && 
							min.pStart.plus <= (res1[k] + this.pattern.size -1)) || 
							 (max.pStart.plus >= res1[k] && 
							  max.pStart.plus <= (res1[k] + this.pattern.size -1))))
                            c(as.character(seqs[i]),seq.names[i], pattern.name,
                                as.character(reverseComplement(patterns[[j]])), 
                                this.pattern.size - 1 + res1[k], res1[k])
						else if (res1[k] >0 && overlap.allpos && 
							min.pStart.plus >= res1[k] && 
							max.pStart.plus <= (res1[k] + this.pattern.size -1))
							c(as.character(seqs[i]),seq.names[i], pattern.name,
								as.character(reverseComplement(patterns[[j]])), 
							    this.pattern.size - 1 + res1[k], res1[k])
                   }))
                }
            ))
            if (length(minus.gRNAs) >1)
            {
                #print(minus.gRNAs)
                # print(j)
                # print(dim(minus.gRNAs))
                colnames(minus.gRNAs) <- c("gRNAPlusPAM", 
                    "REcutgRNAName", "REname", "REpattern",
                    "REcutStart","REcutEnd")
                gRNAs.RE <- rbind(gRNAs.RE,  minus.gRNAs)	
            }
        }		 
        pattern <- translatePattern(pattern)
        pos.gRNAs <- do.call(rbind, lapply(1:length(all.gRNAs), 
            function(i){
                res1  <- as.numeric(gregexpr(pattern, seqs[i], perl = TRUE,
                    fixed = FALSE, ignore.case = TRUE)[[1]])
                do.call(rbind, lapply(1:length(res1), function(k) {
					if (res1[k] >0 && !overlap.allpos && 
						((min.pStart.plus >= res1[k] && 
						min.pStart.plus <= (res1[k] + this.pattern.size -1)) || 
						(max.pStart.plus >= res1[k] && 
						max.pStart.plus <= (res1[k] + this.pattern.size -1))))
                        c(as.character(seqs[i]), seq.names[i], pattern.name, 
                            as.character(patterns[[j]]), res1[k], 
                            res1[k] + this.pattern.size - 1)
					else if (res1[k] >0 && overlap.allpos && 
						min.pStart.plus >= res1[k] && 
						max.pStart.plus <= (res1[k] + this.pattern.size -1))
						c(as.character(seqs[i]), seq.names[i], pattern.name, 
						    as.character(patterns[[j]]), res1[k], 
							res1[k] + this.pattern.size - 1)
                }))
            }
        ))
        if (length(pos.gRNAs) > 0)
        {
            #print("pos.gRNAs")
            #print(j)
            #print(pos.gRNAs)
            colnames(pos.gRNAs) <- c("gRNAPlusPAM", "REcutgRNAName", 
                "REname", "REpattern", "REcutStart", "REcutEnd")
            gRNAs.RE <- rbind(gRNAs.RE,  pos.gRNAs)	
        }
    }
    gRNAs.RE <- gRNAs.RE[gRNAs.RE[,1] != "", ]
    if (! missing(pairOutputFile) && file.exists(pairOutputFile)) {
        gRNAs.RE.plus <- gRNAs.RE
        gRNAs.RE.minus <- gRNAs.RE
        colnames(gRNAs.RE.plus) <- paste("Forward", 
            colnames(gRNAs.RE.plus), sep = "")
        colnames(gRNAs.RE.minus) <- paste("Reverse", 
            colnames(gRNAs.RE.minus), sep = "")
        pairgRNAs <- read.table(pairOutputFile, sep = "\t", header = TRUE,
            stringsAsFactors = FALSE)
        ann.gRNAs <- merge(pairgRNAs, gRNAs.RE.plus, 
            by = "ForwardgRNAPlusPAM", all.x = TRUE)
        ann.gRNAs <- merge(ann.gRNAs, gRNAs.RE.minus, 
            by = "ReversegRNAPlusPAM", all.x = TRUE)
        ann.gRNAs <- cbind(ann.gRNAs[,1], ann.gRNAs[,3], ann.gRNAs[,2],
            ann.gRNAs[,4:dim(ann.gRNAs)[2]])
        colnames(ann.gRNAs)[1:3] <- c("ReversegRNAPlusPAM", 
            "ReversegRNAName", "ForwardgRNAPlusPAM")
        ann.gRNAs <- ann.gRNAs[order(ann.gRNAs[,1]), ]
        withRE <- ann.gRNAs[ ! is.na(ann.gRNAs$ForwardREname) | 
            ! is.na(ann.gRNAs$ReverseREname), ]
        withRE <- unique(cbind(withRE[,1:4]))
        if (dim(withRE)[1] == 0 && findgRNAsWithREcutOnly)
            stop("No pairs with RE sites!")
        gRNAs  <- DNAStringSet(c(as.character(withRE$ForwardgRNAPlusPAM),
            as.character(withRE$ReversegRNAPlusPAM)))
        names(gRNAs) <- c(as.character(withRE$ForwardgRNAName),
            as.character(withRE$ReversegRNAName))
    }
    else ### from unpaired search
    {
        if (dim(gRNAs.RE)[1] == 0 && findgRNAsWithREcutOnly)
            stop("No gRNAs with RE sites!")
        temp <- unique(cbind(as.character(gRNAs.RE$gRNAPlusPAM), 
            as.character(gRNAs.RE$REcutgRNAName)))
        gRNAs <- DNAStringSet(temp[, 1])
        names(gRNAs) <- temp[, 2]
        ann.gRNAs <- gRNAs.RE
    }
    list(gRNAs.withRE = unique(gRNAs), gRNAREcutDetails = ann.gRNAs)
}
