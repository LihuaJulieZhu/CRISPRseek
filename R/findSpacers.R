findSpacers <-
    function (inputFilePath, format = "fasta", PAM = "NGG", PAM.size = 3, 
        findPairedSpacerOnly = FALSE, spacer.size = 20, min.gap = 0, 
        max.gap = 20, pairOutputFile, name.prefix = "gRNA") 
{
    if (missing(inputFilePath)) {
        stop("inputFilePath containing the searching sequence is required!")
    }
    if (! file.exists(inputFilePath)) {
        stop("inputfile specified as ", inputFilePath, " does not exists!")
    }
    if (format != "fasta" && format != "fastq") 
    {
        stop("format needs to be either fasta or fastq!")
    }
    subjects <- readDNAStringSet(inputFilePath, format, use.names = TRUE)
    if (width(subjects) == 0)
    {
        stop("The input file contains no sequence! This could be caused by 
            wrong format of the file. If file is created in mac, you could 
            reformat to text by typing tr \"\\r\" \"\\n\" >newfile in the 
            command line")
    }
    for (i in 1:length(subjects))
    {
        subjectname <- gsub("'", "", names(subjects)[i])
        subjectname <- gsub(" ", "", subjectname)
        subject <- subjects[[i]]
        revsubject <- reverseComplement(subject)
        PAM <- gsub("N", "[A|C|G|T]", PAM)
        PAM <- paste(PAM, "$", sep = "")
        seq.len <- nchar(as.character(subject))
        x <- seq.len - spacer.size - PAM.size + 1
        pos.spacers <- do.call(rbind, lapply(1:x, function(i){
            seq <- subseq(as.character(subject), i, 
                (i + spacer.size + PAM.size -1))
            if (regexpr(PAM, seq, perl = TRUE,ignore.case = TRUE)[[1]] > 0)
                c(seq, paste(subjectname,"Start",i,"End", 
                    (i + spacer.size + PAM.size-1), sep = ""),i,"+")
        }))
        minus.spacers <- do.call(rbind, lapply(1:x, function(i){
            seq <- subseq(as.character(revsubject), i, 
                (i + spacer.size + PAM.size -1))
            if (regexpr(PAM, seq, perl = TRUE,ignore.case = TRUE)[[1]] > 0)
                c(seq, paste( subjectname,"Start", (seq.len - i + 1), "End", 
                    (seq.len - i + 1 - spacer.size - PAM.size + 1), sep = ""),
                    (seq.len - i + 1), "-")			
        }))
        plus.index <- numeric()
        minus.index <- numeric()
        forward.index <- 0
        reverse.index <- 0
        if (length(pos.spacers) > 1)
            n.plus.spacers <- dim(pos.spacers)[1]
        else
            n.plus.spacers <- 0
        if (length(minus.spacers) > 1)
            n.minus.spacers <- dim(minus.spacers)[1]
        else
            n.minus.spacers <- 0
        if (findPairedSpacerOnly && (n.minus.spacers  * n.plus.spacers) == 0)
        {
            stop("No paired spacers found in the input sequence!")
        }
        temp = matrix(nrow =0, ncol=5)
        colnames(temp)[1:5] <- c( "ReverseSpacerPlusPAM",
            "ReverseSpacerName", "ForwardSpacerPlusPAM",
            "ForwardSpacerName", "gap")
        write.table(temp, file = pairOutputFile, sep = "\t",
            row.names = FALSE, quote = FALSE)
        if (n.minus.spacers > 0 && n.plus.spacers > 0)
        {
            for (j in 1:n.plus.spacers)
            {
                for (k in 1:n.minus.spacers)
                {
                    if ((as.numeric(as.character(pos.spacers[j,3])) - 1 - 
                        as.numeric(as.character(minus.spacers[k,3]))) >= min.gap
                        && (as.numeric(as.character(pos.spacers[j,3])) - 
                        as.numeric(as.character(minus.spacers[k,3]))) <= 
                        max.gap)
                    {
                        if(findPairedSpacerOnly)
                        {
                            if (!k %in% minus.index)
                            {
                                reverse.index <- reverse.index + 1
                                minus.spacers[k,2] <- paste(name.prefix, "r", 
                                    reverse.index, "_", minus.spacers[k, 2], 
                                    sep = "")
                            }
                            if (!j %in% plus.index)
                            {
                                forward.index <- forward.index + 1
                                pos.spacers[j,2] <- paste(name.prefix,"f",
                                    forward.index, "_", pos.spacers[j, 2], 
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
                paired <- cbind(minus.spacers[minus.index,1], 
                    minus.spacers[minus.index, 2], pos.spacers[plus.index,1], 
                    pos.spacers[plus.index, 2],
                    gap = as.numeric(as.character(
                    pos.spacers[plus.index,3])) - 1 - as.numeric(
                    as.character(minus.spacers[minus.index, 3])))
                colnames(paired)[1:5] <- c( "ReverseSpacerPlusPAM", 
                    "ReverseSpacerName", "ForwardSpacerPlusPAM", 
                    "ForwardSpacerName", "gap")
                if(findPairedSpacerOnly)
                {
                    plus.index <- unique(plus.index)
                    minus.index <- unique(minus.index)
                    all.spacers <- DNAStringSet(c(pos.spacers[plus.index, 1], 
                        minus.spacers[minus.index,1]))
                    names(all.spacers) <- c(pos.spacers[plus.index,2], 
                        minus.spacers[minus.index,2])
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
                        row.names = FALSE, quote = FALSE)
                else
                    write.table(paired[order(as.character(paired[,4])), ], 
                        file = pairOutputFile, sep = "\t", row.names = FALSE,
                        quote = FALSE)
            }
            else if (findPairedSpacerOnly)
            {
                stop("No paired spacers found!")
            }
        }### if pos.spacers and minus.spacers not empty
        if (! findPairedSpacerOnly)
        {
            all.spacers <- DNAStringSet(c(pos.spacers[,1], minus.spacers[,1]))
            if (length(all.spacers) == 0)
                stop("No spacers found in the input sequence!")
            if (n.plus.spacers == 0)
                names(all.spacers) <- paste(name.prefix,
                    c(rep("r", dim(minus.spacers)[1])),
                    c(1:dim(minus.spacers)[1]), "_", 
                    c(minus.spacers[,2]), sep = "")
            else if (n.minus.spacers == 0)
                names(all.spacers) <- paste(name.prefix,c(rep("f", 
                    dim(pos.spacers)[1])),
                    c(1:dim(pos.spacers)[1]), 
                    c(pos.spacers[,2]), sep = "")
            else
                names(all.spacers) <- paste(name.prefix,c(rep("f", 
                    dim(pos.spacers)[1]), rep("r", dim(minus.spacers)[1])),
                    c(1:dim(pos.spacers)[1], 1:dim(minus.spacers)[1]), "_", 
                    c(pos.spacers[,2], minus.spacers[,2]), sep = "")
        }
    } ## for subjects
    all.spacers
}
