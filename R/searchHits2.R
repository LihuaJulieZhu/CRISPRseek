.preprocess_me2 <- function(gRNAs, max.mismatch)
{
    if (length(gRNAs) <= 10L)
        return(FALSE)
    gRNA_min_width <- min(width(gRNAs))
    if (gRNA_min_width < 6L || gRNA_min_width %/% (max.mismatch + 1L) < 3L)
        return(FALSE)
    if (sum(alphabetFrequency(gRNAs, baseOnly=TRUE)[ , "other"]) != 0L)
        return(FALSE)
    TRUE
}

.searchHitsInOneSeq2 <- function(gRNAs, seq, seqname, PAM, PAM.pattern, PAM.size,
                                max.mismatch, outfile, allowed.mismatch.PAM,
                                PAM.location = "3prime", BSgenomeName,
                                baseEditing = FALSE, targetBase = "C", editingWindow = 4:8)
{
    if (.preprocess_me2(gRNAs, max.mismatch)) {
        patterns <- PDict(gRNAs, max.mismatch = max.mismatch)
    } else {
        patterns <- gRNAs
    }
    all_plus_matches <- matchPDict(patterns, seq,
        max.mismatch = max.mismatch)
    revseq <- reverseComplement(seq)
    all_minus_macthes <- matchPDict(patterns, revseq,
        max.mismatch = max.mismatch)
    for (i in seq_len(length(gRNAs))) {
        patternID <- gsub("'", "", names(gRNAs)[i])
        if (length(patternID) < 1) {
            patternID <- paste("pattern", i, sep = "")
        }
        pattern <- gRNAs[[i]]
        ### by default PAM is NGG or NAG
        plus_matches <- Views(seq, all_plus_matches[[i]])
        if (length(plus_matches) > 0) {
            names(plus_matches) <- rep.int(patternID, length(plus_matches))
            writeHits2(gRNA = pattern, seqname = seqname, 
                matches = plus_matches, strand = "+",
                file = outfile, gRNA.size = length(pattern),
                PAM = PAM, PAM.pattern = PAM.pattern, 
                max.mismatch = max.mismatch,
                chrom.len = length(seq), append = TRUE, 
                PAM.location = PAM.location,
                PAM.size = PAM.size,
                allowed.mismatch.PAM = allowed.mismatch.PAM,
                BSgenomeName = BSgenomeName,
                baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow)
        }
        if (reverseComplement(pattern) != pattern) {
            minus_matches <- Views(revseq, all_minus_macthes[[i]])
            if (length(minus_matches) > 0) {
                names(minus_matches) <- rep.int(patternID,
                    length(minus_matches))
                writeHits2(gRNA = pattern, seqname = seqname,
                    matches = minus_matches, strand = "-",
                    file = outfile, gRNA.size = length(pattern),
                    PAM = PAM, PAM.pattern = PAM.pattern, 
                    max.mismatch = max.mismatch,
                    chrom.len = length(seq), append = TRUE, 
                    PAM.location = PAM.location, PAM.size = PAM.size,
                    allowed.mismatch.PAM = allowed.mismatch.PAM,
                    BSgenomeName = BSgenomeName,
                    baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow) 
            }
        }
    }
}

searchHits2 <-
    function (gRNAs, BSgenomeName, chromToSearch = "all", chromToExclude = "", 
	max.mismatch = 3, PAM.size = 3, gRNA.size = 20, 
        PAM = "NGG", PAM.pattern = "N[A|G]G$",
        allowed.mismatch.PAM = 1, PAM.location = "3prime", 
        baseEditing = FALSE, targetBase = "C", editingWindow = 4:8)  
{
    if (missing(gRNAs) || class(gRNAs) != "DNAStringSet") {
        stop("gRNAs is required as a DNAStringSet object!")
    }
    if (missing(BSgenomeName) || class(BSgenomeName) != "BSgenome") {
        stop("BSgenomeName is required as BSgenome object!")
    }
    outfile <- tempfile(tmpdir = getwd())
    max.mismatch <- max.mismatch
    seqnames <- seqnames(BSgenomeName)
    if (chromToSearch != "all")
        seqnames <- intersect(seqnames, chromToSearch)
    if (length(chromToExclude) >0)
	seqnames <- setdiff(seqnames, chromToExclude)
    append <- FALSE
    if (width(gRNAs)[1] == (gRNA.size + PAM.size))
    {
        if(PAM.location == "3prime")
            gRNAs <- subseq(gRNAs, 1, gRNA.size)
        else
            gRNAs <- subseq(gRNAs, PAM.size + 1, gRNA.size + PAM.size)
    }
    else if (width(gRNAs)[1] != gRNA.size)
    {
        stop("the gRNA length needs to be equal to the 
            specified gRNA.size (or gRNA.size plus PAM.size\n)");
    }
    for (seqname in seqnames) {
        cat(">>> Finding all hits in sequence", seqname, "...\n")
        subject <- BSgenomeName[[seqname]]
        .searchHitsInOneSeq2(gRNAs = gRNAs, seq = subject, 
             seqname = seqname, PAM = PAM, 
             PAM.pattern = PAM.pattern, PAM.size = PAM.size,
             max.mismatch = max.mismatch, outfile, 
             allowed.mismatch.PAM = allowed.mismatch.PAM,
             PAM.location = PAM.location,
             BSgenomeName = BSgenomeName,
             baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow)
        cat(">>> DONE searching\n")
    }
    if (file.exists(outfile))
    {
        hits <- read.table(outfile, sep="\t", header = TRUE, 
            stringsAsFactors = FALSE)
        unlink(outfile)
        hits
    }
    else
    {
        warning("No matching found, please check your input sequence, and make
            sure you are using the right genome. You can also alter your 
            search criteria such as increasing max.mismatch!")
        data.frame()
    }
}
