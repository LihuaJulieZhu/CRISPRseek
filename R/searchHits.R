searchHits <-
    function (gRNAs, BSgenomeName, chromToSearch = "all", max.mismatch = 4, 
        PAM.size = 3, gRNA.size = 20, PAM = "N[A|G]G$") 
{
    if (missing(gRNAs) || class(gRNAs) != "DNAStringSet") {
        stop("dictO is required as a DNAStringSet object!")
    }
    if (missing(BSgenomeName) || class(BSgenomeName) != "BSgenome") {
        stop("BSgenomeName is required as BSgenome object!")
    }
    outfile <- tempfile(tmpdir = getwd())
    max.mismatch <- max.mismatch + PAM.size -1
    seqnames <- seqnames(BSgenomeName)
    if (chromToSearch != "all")
        seqnames <- intersect(seqnames, chromToSearch)
    append <- FALSE
    for (seqname in seqnames) {
        subject <- BSgenomeName[[seqname]]
        revsubject <- reverseComplement(subject)
        chrom.len <- length(subject)
        cat(">>> Finding all hits in sequences", seqname, "...\n")
        for (i in seq_len(length(gRNAs))) {
            patternID <- gsub("'", "", names(gRNAs)[i])
            if (length(patternID) < 1) {
                patternID <- paste("pattern", i, sep = "")
            }
            pattern <- DNAString(toupper(gRNAs[[i]]))
            ### by default PAM is NGG or NAG
            plus_matches <- matchPattern(pattern, subject,
                max.mismatch = max.mismatch, min.mismatch = 0, 
                with.indels = FALSE, fixed = TRUE, algorithm = "auto")
            if (length(plus_matches) > 0) {
                names(plus_matches) <- rep.int(patternID, length(plus_matches))
                writeHits(pattern, seqname, plus_matches, strand = "+", 
                    file = outfile, gRNA.size = gRNA.size,
                    PAM = PAM, max.mismatch = max.mismatch - 2,
                    chrom.len = chrom.len, append = append)
                append <- TRUE
            }            
            if (reverseComplement(pattern) != pattern) {
                minus_matches <- matchPattern(pattern, revsubject, 
                    max.mismatch = max.mismatch, min.mismatch = 0, 
                    with.indels = FALSE, fixed = TRUE, algorithm = "auto")
                if (length(minus_matches) > 0) {
                    names(minus_matches) <- rep.int(patternID,
                        length(minus_matches))
                    writeHits(pattern, seqname, minus_matches, strand = "-", 
                        file = outfile, gRNA.size = gRNA.size, PAM = PAM,
                        max.mismatch = max.mismatch - 2, 
                        chrom.len = chrom.len, append = append)
                    append <- TRUE
                }
            }
        }
        cat(">>> DONE searching\n")
    }
    if (file.exists(outfile))
    {
        hits <- read.table(outfile, sep="\t", header = TRUE, stringsAsFactors = FALSE)
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
