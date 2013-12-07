filterSpacers <-
    function(all.spacers, pairOutputFile = "", REpatternFile, format = "fasta",
        minREpatternSize = 6, overlap.spacer.positions = c(17, 18))
{
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
    seqs <- as.character(all.spacers)
    seq.names <- names(all.spacers)
    min.pStart.plus <- overlap.spacer.positions[1] - width(patterns) + 1
    max.pStart.plus <- overlap.spacer.positions[2] + width(patterns) - 1
    spacers.RE <- data.frame(SpacerPlusPAM = "", REcutSpacerName = "", 
        REname = "", REpattern = "", REcutStart = "", REcutEnd = "")
    for (j in 1:length(patterns))
    {
        pattern.name <- gsub("'", "", names(patterns)[j])
        pattern <- patterns[[j]]
        revpattern <- reverseComplement(pattern)
        if (revpattern != pattern)
        {
            revpattern <- translatePattern(revpattern)
            minus.spacers <- do.call(rbind, lapply(1:length(all.spacers), 
                function(i){
                    res1  <- as.numeric(gregexpr(revpattern, seqs[i],
                        perl = TRUE,fixed = FALSE,ignore.case = TRUE)[[1]])
                    do.call(rbind, lapply(1:length(res1), function(k) {
                        if (res1[k] >0 && res1[k] >= min.pStart.plus[j] && 
                            res1[k] <= max.pStart.plus[j])
                            c(as.character(seqs[i]),seq.names[i], pattern.name,
                                as.character(reverseComplement(patterns[[j]])), 
                                length(patterns[[j]]) - 1 + res1[k], res1[k])			
                   }))
                }
            ))
            if (length(minus.spacers) >1)
            {
                #print(minus.spacers)
                # print(j)
                # print(dim(minus.spacers))
                colnames(minus.spacers) <- c("SpacerPlusPAM", 
                    "REcutSpacerName", "REname", "REpattern",
                    "REcutStart","REcutEnd")
                spacers.RE <- rbind(spacers.RE,  minus.spacers)	
            }
        }		 
        pattern <- translatePattern(pattern)
        pos.spacers <- do.call(rbind, lapply(1:length(all.spacers), 
            function(i){
                res1  <- as.numeric(gregexpr(pattern, seqs[i], perl = TRUE,
                    fixed = FALSE, ignore.case = TRUE)[[1]])
                do.call(rbind, lapply(1:length(res1), function(k) {
                    if (res1[k] >0 && res1[k] >= min.pStart.plus[j] && 
                        res1[k] <= max.pStart.plus[j])
                        c(as.character(seqs[i]), seq.names[i], pattern.name, 
                            as.character(patterns[[j]]), res1[k], 
                            res1[k] + length(patterns[[j]])- 1)			
                }))
            }
        ))
        if (length(pos.spacers) > 0)
        {
            #print("pos.spacers")
            #print(j)
            #print(pos.spacers)
            colnames(pos.spacers) <- c("SpacerPlusPAM", "REcutSpacerName", 
                "REname", "REpattern", "REcutStart", "REcutEnd")
            spacers.RE <- rbind(spacers.RE,  pos.spacers)	
        }
    }
    spacers.RE <- spacers.RE[spacers.RE[,1] != "", ]
    if (! missing(pairOutputFile) && file.exists(pairOutputFile)) {
        spacers.RE.plus <- spacers.RE
        spacers.RE.minus <- spacers.RE
        colnames(spacers.RE.plus) <- paste("Forward", 
            colnames(spacers.RE.plus), sep = "")
        colnames(spacers.RE.minus) <- paste("Reverse", 
            colnames(spacers.RE.minus), sep = "")
        pairSpacers <- read.table(pairOutputFile, sep = "\t", header = TRUE,
            stringsAsFactors = FALSE)
        ann.Spacers <- merge(pairSpacers, spacers.RE.plus, 
            by = "ForwardSpacerPlusPAM", all.x = TRUE)
        ann.Spacers <- merge(ann.Spacers, spacers.RE.minus, 
            by = "ReverseSpacerPlusPAM", all.x = TRUE)
        ann.Spacers <- cbind(ann.Spacers[,1], ann.Spacers[,3], ann.Spacers[,2],
            ann.Spacers[,4:dim(ann.Spacers)[2]])
        colnames(ann.Spacers)[1:3] <- c("ReverseSpacerPlusPAM", 
            "ReverseSpacerName", "ForwardSpacerPlusPAM")
        ann.Spacers <- ann.Spacers[order(ann.Spacers[,1]), ]
        withRE <- ann.Spacers[ ! is.na(ann.Spacers$ForwardREname) | 
            ! is.na(ann.Spacers$ReverseREname), ]
        withRE <- unique(cbind(withRE[,1:4]))
        if (dim(withRE)[1] == 0)
            stop("No pairs with RE sites!")
        spacers  <- DNAStringSet(c(as.character(withRE$ForwardSpacerPlusPAM),
            as.character(withRE$ReverseSpacerPlusPAM)))
        names(spacers) <- c(as.character(withRE$ForwardSpacerName),
            as.character(withRE$ReverseSpacerName))
    }
    else ### from unpaired search
    {
        if (dim(spacers.RE)[1] == 0)
            stop("No spacers with RE sites!")
        temp <- unique(cbind(as.character(spacers.RE$SpacerPlusPAM), 
            as.character(spacers.RE$REcutSpacerName)))
        spacers <- DNAStringSet(temp[, 1])
        names(spacers) <- temp[, 2]
        ann.Spacers <- spacers.RE
    }
    list(spacers.withRE = unique(spacers), spacerREcutDetails = ann.Spacers)
}
