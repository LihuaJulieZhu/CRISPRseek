.preprocess_me <- function(gRNAs, max.mismatch)
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

.searchHitsInOneSeq <- function(gRNAs, seq, seqname, PAM, PAM.size,
                                max.mismatch, outfile)
{
    if (.preprocess_me(gRNAs, max.mismatch)) {
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
            writeHits(pattern, seqname, plus_matches, strand = "+",
                file = outfile, gRNA.size = length(pattern) - PAM.size,
                PAM = PAM, max.mismatch = max.mismatch - 2,
                chrom.len = length(seq), append = TRUE)
        }
        if (reverseComplement(pattern) != pattern) {
            minus_matches <- Views(revseq, all_minus_macthes[[i]])
            if (length(minus_matches) > 0) {
                names(minus_matches) <- rep.int(patternID,
                    length(minus_matches))
                writeHits(pattern, seqname, minus_matches, strand = "-",
                    file = outfile, gRNA.size = length(pattern) - PAM.size,
                    PAM = PAM, max.mismatch = max.mismatch - 2,
                    chrom.len = length(seq), append = TRUE)
            }
        }
    }
}

searchHits <-
    function (gRNAs, BSgenomeName, chromToSearch = "all", max.mismatch = 4, 
        PAM.size = 3, gRNA.size = 20, PAM = "N[A|G]G$") 
{
    if (missing(gRNAs) || class(gRNAs) != "DNAStringSet") {
        stop("gRNAs is required as a DNAStringSet object!")
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
        cat(">>> Finding all hits in sequence", seqname, "...\n")
        subject <- BSgenomeName[[seqname]]
        .searchHitsInOneSeq(gRNAs, subject, seqname, PAM, PAM.size,
                            max.mismatch, outfile)
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
