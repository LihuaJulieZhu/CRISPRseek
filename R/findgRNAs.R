.getgRNA.cut.sites <- function(subject, subjectname, PAM ="NGG", 
    gRNA.pattern = "", gRNA.size = 20, cut.site = 17,
    PAM.size = 3, calculategRNAEfficacy = TRUE, baseBeforegRNA = 4,
    baseAfterPAM = 3,
    reverse.subject = FALSE,
    PAM.location = "3prime")
{
     seq.len <- nchar(as.character(subject))
     pos.PAMs <- unlist(gregexpr(PAM, subject, perl = TRUE,
         ignore.case = TRUE, fixed = FALSE))
     if (PAM.location == "3prime")
         pos.PAMs <- pos.PAMs[pos.PAMs != -1 & pos.PAMs > gRNA.size]
     else
         pos.PAMs <- pos.PAMs[pos.PAMs != -1]
     if (length(pos.PAMs) > 0)
     {
        if (PAM.location == "3prime")
            starts.gRNA <- pos.PAMs - gRNA.size
        else
            starts.gRNA <- pos.PAMs + PAM.size 
        starts.gRNA <- subset(starts.gRNA, starts.gRNA <= (seq.len - gRNA.size + 1)) 
        ends.gRNA <- starts.gRNA + gRNA.size - 1
        if (gRNA.pattern != "")
        {
            gRNA.seqs <- as.character(Views(subject, 
               start = starts.gRNA,
               end = ends.gRNA))
            pos.set2 <- unlist(gregexpr(gRNA.pattern, gRNA.seqs,
               perl = TRUE, ignore.case = TRUE, fixed = FALSE))
            pos.PAMs <- pos.PAMs[pos.set2 == 1]
        }
        if (PAM.location == "3prime")
            seq <- as.character(Views(subject,
                start = starts.gRNA,
                end = ends.gRNA + PAM.size))
        else
            seq <- as.character(Views(subject,
                start = starts.gRNA - PAM.size,
                end = ends.gRNA))
        extendedSequence <- seq
        if (calculategRNAEfficacy)
        {
            extended.starts <- starts.gRNA - baseBeforegRNA
            extended.starts[extended.starts < 0] <- 1
            if (PAM.location == "3prime")
                extended.ends <- ends.gRNA + PAM.size + baseAfterPAM
            else
                extended.ends <- starts.gRNA + baseAfterPAM - 1
            extended.ends[extended.ends > length(subject)] <- length(subject)
            extendedSequence <- as.character(Views(subject, start = extended.starts,
                end = extended.ends))
        }
        if (reverse.subject)
        {
            gRNAs.cut <- cbind(seq, paste( subjectname,"_gR",
                (seq.len - (starts.gRNA + cut.site -1 ) + 1), "r", sep = ""),
                (seq.len - starts.gRNA + 1), "-", extendedSequence)
        }
        else
        {
            gRNAs.cut <-
                cbind(seq, paste(subjectname,"_gR", (starts.gRNA + cut.site - 1),
                "f", sep = ""), starts.gRNA, "+", extendedSequence)
        }
     }
     else
     {
        gRNAs.cut <- ""
     }
        gRNAs.cut
}

## Very inefficient (perform in quadratic time).
.compute_pair_index <- function(plus_start, minus_start, min.gap, max.gap)
{
    plus.index <- minus.index <- integer(0)
    for (j in seq_along(plus_start))
    {
        for (k in seq_along(minus_start))
        {
            d <- plus_start[j] - minus_start[k]
            if (d > min.gap && d <= max.gap)
            {
                plus.index <- c(plus.index, j)
                minus.index <- c(minus.index, k)
            } ### if paired
        } #for k
    }## for j
    list(plus.index, minus.index)
}

## A more efficient version of .compute_pair_index() that performs almost in
## linear time.
.compute_pair_index2 <- function(plus_start, minus_start, min.gap, max.gap)
{
    plus.index <- minus.index <- integer(0)
    n1 <- length(plus_start)
    n2 <- length(minus_start)
    if (n1 == 0L || n2 == 0L)
        return(list(plus.index, minus.index))
    plus_oo <- order(plus_start)
    minus_oo <- order(minus_start)
    i2 <- 1L
    for (i1 in seq_len(n1)) {
        j <- plus_oo[[i1]]
        current_plus_start <- plus_start[[j]]
        min_minus_start <- current_plus_start - max.gap
        max_minus_start <- current_plus_start - min.gap - 1
        while (minus_start[[minus_oo[[i2]]]] < min_minus_start) {
            i2 <- i2 + 1L
            if (i2 > n2)
                break
        }
        if (i2 > n2)
            break
        i <- i2
        while (minus_start[[k <- minus_oo[[i]]]] <= max_minus_start) {
            plus.index <- c(plus.index, j)
            minus.index <- c(minus.index, k)
            i <- i + 1L
            if (i > n2)
                break
        }
    }
    oo <- orderIntegerPairs(plus.index, minus.index)
    list(plus.index[oo], minus.index[oo])
}

## An even more efficient version of .compute_pair_index() based on
## findOverlaps(). Hundreds times faster than .compute_pair_index2().
.compute_pair_index3 <- function(plus_start, minus_start, min.gap, max.gap)
{
    shift <- floor((min.gap + max.gap) / 2)
    subject <- IRanges(minus_start + shift, width=1L)
    maxgap <- ceiling((max.gap - min.gap) / 2)
    query <- IRanges(plus_start, width=1L)
    hits <- findOverlaps(query, subject, maxgap=maxgap)
    d <- plus_start[queryHits(hits)] - minus_start[subjectHits(hits)]
    hits <- sort(hits[min.gap < d & d <= max.gap])
    list(queryHits(hits), subjectHits(hits))
}

findgRNAs <-
    function (inputFilePath, format = "fasta", PAM = "NGG", PAM.size = 3,
        findPairedgRNAOnly = FALSE, annotatePaired = TRUE,
        enable.multicore = FALSE, n.cores.max = 6,
        gRNA.pattern = "", gRNA.size = 20, 
	overlap.gRNA.positions = c(17,18),
        min.gap = 0, max.gap = 20, pairOutputFile, name.prefix = "",
	    featureWeightMatrixFile = system.file("extdata", 
        "DoenchNBT2014.csv", package = "CRISPRseek"), baseBeforegRNA = 4, 
	    baseAfterPAM = 3,
	    calculategRNAEfficacy = FALSE, efficacyFile,
         PAM.location = "3prime") 
{
    if (missing(inputFilePath)) {
        stop("inputFilePath containing the searching sequence 
	   or a DNAStringSet object is required!")
    }
    cut.site <- min(overlap.gRNA.positions)
    if (class(inputFilePath) != "DNAStringSet")
    {
    	if (! file.exists(inputFilePath)) {
           stop("inputfile specified as ", inputFilePath, " does not exists!")
        }
        if (format == "fasta" || format == "fastq")
        {
            subjects <- readDNAStringSet(inputFilePath, format, use.names = TRUE)
        }
        else
        {
            stop("format needs to be either fasta or fastq !")
        }
    }
    else
    {
	    subjects <- inputFilePath
    }
    PAM <- translatePattern(PAM)
    PAM <- paste("(?=", PAM, ")", sep="")

    gRNA.pattern <- translatePattern(gRNA.pattern)
    min.subject <- gRNA.size + PAM.size 
    subjects <-  subjects[width(subjects) >= min.subject,]
    if (length(subjects) == 0)
    {
        stop("The input file contains no sequence! This could be caused by 
            wrong format of the file. If file is created in mac, you could 
            reformat to text by typing tr \"\\r\" \"\\n\" >newfile in the 
            command line")
    }
    toAppend = FALSE
    colNames = TRUE
    names(subjects) <- gsub("'", "", names(subjects))
    names(subjects) <- gsub(" ", "", names(subjects))
    names(subjects) <- gsub("\t", ":", names(subjects))
     
    all.gRNAs.df <- do.call(rbind, lapply(1:length(subjects), function(p)
   {
        subject <- subjects[[p]]
        subjectname <- names(subjects)[p]
        revsubject <- reverseComplement(subject)
        plus.gRNAs <- 
           .getgRNA.cut.sites(subject, subjectname, PAM = PAM, 
               gRNA.pattern = gRNA.pattern, 
               gRNA.size = gRNA.size,
               cut.site = cut.site,
               PAM.size = PAM.size, 
               calculategRNAEfficacy = calculategRNAEfficacy,
               baseBeforegRNA = baseBeforegRNA,
               baseAfterPAM = baseAfterPAM,
               PAM.location = PAM.location)
       
       minus.gRNAs <-
           .getgRNA.cut.sites(revsubject, subjectname, PAM = PAM,
               gRNA.pattern = gRNA.pattern,
               gRNA.size = gRNA.size,
               cut.site = cut.site,
               PAM.size = PAM.size,
               calculategRNAEfficacy = calculategRNAEfficacy,
               baseBeforegRNA = baseBeforegRNA,
               baseAfterPAM = baseAfterPAM,
               reverse.subject = TRUE,
               PAM.location = PAM.location) 

        if (length(plus.gRNAs) > 1)
            n.plus.gRNAs <- dim(plus.gRNAs)[1]
        else
            n.plus.gRNAs <- 0
        if (length(minus.gRNAs) > 1)
            n.minus.gRNAs <- dim(minus.gRNAs)[1]
        else
            n.minus.gRNAs <- 0
        if (n.minus.gRNAs == 0 && n.plus.gRNAs == 0)
        {
            cat(paste("No gRNAs found in the input sequence", subjectname))
            all.gRNAs <- DNAStringSet()
        }
        else if (findPairedgRNAOnly && (n.minus.gRNAs  * n.plus.gRNAs) == 0)
        {
            cat(paste("No paired gRNAs found in the input sequence", subjectname))
            all.gRNAs <- DNAStringSet()
        }
	    else if (annotatePaired || findPairedgRNAOnly)
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
                plus_start <- as.numeric(as.character(plus.gRNAs[ , 3L]))
                minus_start <- as.numeric(as.character(minus.gRNAs[ , 3L]))
                pair_index <- .compute_pair_index3(plus_start,
                                                   minus_start,
                                                   min.gap, max.gap)
                plus.index <- pair_index[[1L]]
                minus.index <- pair_index[[2L]]
                if (!findPairedgRNAOnly)
                {
                    all.gRNAs <- DNAStringSet(c(plus.gRNAs[,1], minus.gRNAs[,1]))
                    names(all.gRNAs) <- c(plus.gRNAs[,2], minus.gRNAs[,2])
                    forEffi <- rbind(plus.gRNAs, minus.gRNAs)
                }
                if (length(minus.index) > 0 && length(plus.index) > 0)
                {
                    paired <- cbind(minus.gRNAs[minus.index,1], 
                        minus.gRNAs[minus.index, 2], plus.gRNAs[plus.index,1], 
                        plus.gRNAs[plus.index, 2],
                        gap = as.numeric(as.character(
                        plus.gRNAs[plus.index,3])) - 1 - as.numeric(
                        as.character(minus.gRNAs[minus.index, 3])))
                    colnames(paired)[1:5] <- c( "ReversegRNAPlusPAM", 
                        "ReversegRNAName", "ForwardgRNAPlusPAM", 
                        "ForwardgRNAName", "gap")
                    if(findPairedgRNAOnly)
                    {
                        plus.index <- unique(plus.index)
                        minus.index <- unique(minus.index)
                        all.gRNAs <- DNAStringSet(c(plus.gRNAs[plus.index, 1], 
                            minus.gRNAs[minus.index,1]))
                        names(all.gRNAs) <- c(plus.gRNAs[plus.index,2], 
                            minus.gRNAs[minus.index,2])
                        forEffi <- rbind(plus.gRNAs[plus.index,], minus.gRNAs[minus.index,])
                    }
                    if (dim(paired)[1] == 1)
                        write.table(paired, file = pairOutputFile, sep = "\t", 
                            row.names = FALSE, quote = FALSE, append = toAppend,
			                col.names = colNames)
                    else
                        write.table(paired[order(as.character(paired[,4])), ], 
                            file = pairOutputFile, sep = "\t", row.names = FALSE,
                            quote = FALSE, append = toAppend, col.names = colNames)
                } ## if paired found
                else if (findPairedgRNAOnly)
                {
                    cat(paste("No paired gRNAs found for sequence",
                        subjectname))
                   all.gRNAs <- DNAStringSet()
                   forEffi <- ""
                }
            }### if plus.gRNAs and minus.gRNAs not empty
            else if (!findPairedgRNAOnly && n.plus.gRNAs > 0)
            {
                all.gRNAs <- DNAStringSet(plus.gRNAs[,1])
                names(all.gRNAs) <- plus.gRNAs[,2]
                forEffi <- plus.gRNAs
            }
            else if (!findPairedgRNAOnly && n.minus.gRNAs > 0)
            {
                all.gRNAs <- DNAStringSet(minus.gRNAs[,1])
                names(all.gRNAs) <- minus.gRNAs[,2]
                forEffi <- minus.gRNAs
            }
        } ### annotatePaired and (no paired found or findPairedOnly)
        else if (n.minus.gRNAs > 0 && n.plus.gRNAs > 0)
        {
             all.gRNAs <- DNAStringSet(c(plus.gRNAs[,1], minus.gRNAs[,1]))
             names(all.gRNAs) <- c(plus.gRNAs[,2], minus.gRNAs[,2])
             forEffi <- rbind(plus.gRNAs, minus.gRNAs)
        }
        else if (n.minus.gRNAs > 0 )
        {
             all.gRNAs <- DNAStringSet(minus.gRNAs[,1])
             names(all.gRNAs) <- minus.gRNAs[,2]
             forEffi <- minus.gRNAs
        }
        else if (n.plus.gRNAs > 0 )
        {
            all.gRNAs <- DNAStringSet(plus.gRNAs[,1])
            names(all.gRNAs) <- plus.gRNAs[,2]
            forEffi <- plus.gRNAs
        }
        else
        {
            all.gRNAs <- DNAStringSet()
        }
        if (length(all.gRNAs) == 0 && !findPairedgRNAOnly)
          cat(paste("No gRNAs found in the input sequence",
               subjectname))
        if (length(all.gRNAs) >0)
	       forEffi
    })) ## do call subjects
    if (calculategRNAEfficacy && length(all.gRNAs.df) > 4)
    {
        featureWeightMatrix <- read.csv(featureWeightMatrixFile, header=TRUE)
        effi <- calculategRNAEfficiency(all.gRNAs.df[,5], 
	    baseBeforegRNA = baseBeforegRNA, 
	    featureWeightMatrix = featureWeightMatrix, 
            enable.multicore = enable.multicore,
            n.cores.max = n.cores.max,
            gRNA.size = gRNA.size) 
        extendedSequences <- cbind(all.gRNAs.df, effi)
        colnames(extendedSequences)  <- c("gRNAplusPAM", "name", "start", "strand", 
	    "extendedSequence", "gRNAefficacy")
        if (PAM.location == "3prime")
            extendedSequences[nchar(extendedSequences[,5]) 
	        < baseBeforegRNA + gRNA.size + PAM.size + baseAfterPAM, 6] <- 
                "extended sequence too short"
        else
            extendedSequences[nchar(extendedSequences[,5])
                < baseBeforegRNA + baseAfterPAM , 6] <-
                "extended sequence too short"
        write.table(extendedSequences,
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
