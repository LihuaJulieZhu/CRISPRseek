####### for compare2Sequences
.searchHitsInOneSeq3 <- function(gRNAs, seqs, seqname, 
    PAM = "NGG", PAM.pattern = "NNN$", PAM.size = 3,
    max.mismatch = 2, outfile, allowed.mismatch.PAM = 2,
    PAM.location = "3prime", gRNA.size = 20)
{
    if (.preprocess_me2(gRNAs, max.mismatch)) {
        patterns <- PDict(gRNAs, max.mismatch = max.mismatch)
    } else {
        patterns <- gRNAs
    }
    all_plus_matches <- matchPDict(patterns, seqs,
        max.mismatch = max.mismatch)
    revseq <- reverseComplement(seqs)
    all_minus_matches <- matchPDict(patterns, revseq,
        max.mismatch = max.mismatch)
    for (i in seq_len(length(gRNAs))) {
        patternID <- gsub("'", "", names(gRNAs)[i])
        if (length(patternID) < 1) {
            patternID <- paste("pattern", i, sep = "")
        }
        pattern1 <- gRNAs[[i]]
        ### by default PAM is NGG or NAG
        plus_matches <- Views(seqs, all_plus_matches[[i]])
        if (length(plus_matches) > 0) {
            names(plus_matches) <- rep.int(patternID, length(plus_matches))
            writeHits(gRNA = pattern1, seqname = seqname, 
                matches = plus_matches, strand = "+",
                file = outfile, gRNA.size = gRNA.size,
                PAM = PAM, PAM.pattern = PAM.pattern, 
                max.mismatch = max.mismatch,
                chrom.len = length(seqs), append = TRUE, 
                PAM.location = PAM.location,
                PAM.size = PAM.size,
                allowed.mismatch.PAM = allowed.mismatch.PAM,
                seqs = seqs)
        }
        if (reverseComplement(pattern1) != pattern1) {
            minus_matches <- Views(revseq, all_minus_matches[[i]])
            if (length(minus_matches) > 0) {
                names(minus_matches) <- rep.int(patternID,
                    length(minus_matches))
                writeHits(gRNA = pattern1, seqname = seqname,
                    matches = minus_matches, strand = "-",
                    file = outfile, gRNA.size = gRNA.size,
                    PAM = PAM, PAM.pattern = PAM.pattern, 
                    max.mismatch = max.mismatch,
                    chrom.len = length(seqs), append = TRUE, 
                    PAM.location = PAM.location, PAM.size = PAM.size,
                    allowed.mismatch.PAM = allowed.mismatch.PAM,
                    seqs = revseq)
            }
        }
    }
}

searchHits <-
    function (gRNAs,seqs,seqname,
	max.mismatch = 3, PAM.size = 3, gRNA.size = 20, 
        PAM = "NGG", PAM.pattern = "NNN$",
        allowed.mismatch.PAM = 2, PAM.location = "3prime",
        outfile) 
{
    if (missing(gRNAs) || class(gRNAs) != "DNAStringSet") {
        stop("gRNAs is required as a DNAStringSet object!")
    }
    if (missing(seqs) || class(seqs) != "DNAString") {
        stop("seq is required as a DNAString object!")
    }
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
    cat(">>> Finding all hits in sequence", seqname, "...\n")
        .searchHitsInOneSeq3(gRNAs = gRNAs, seqs = seqs, 
             seqname = seqname,
             PAM = PAM, gRNA.size = gRNA.size,
             PAM.pattern = PAM.pattern, PAM.size = PAM.size,
             max.mismatch = max.mismatch, outfile, 
             allowed.mismatch.PAM = allowed.mismatch.PAM,
             PAM.location = PAM.location)
    cat(">>> DONE searching\n")
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
