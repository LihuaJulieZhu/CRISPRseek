writeHits <-
    function (gRNA, seqname, matches, strand, file, gRNA.size = 20, 
        PAM = "N[A|G]G$", max.mismatch = 4, chrom.len, append = FALSE) 
{
    if (missing(gRNA) || class(gRNA) != "DNAString") {
        stop("gRNA is required as a DNAString object!")
    }	
    if (missing(seqname)) {
        stop("seqname is required as character!")
    }
    if (missing(matches) || class(matches) != "XStringViews") {
        stop("matches is required as XStringViews object!")
    }
    if (missing(strand)) {
        stop("strand is required as + or - !")
    }
    if (file.exists(file) && !append) 
        warning("existing file ", file, 
            " will be overwritten with 'append = FALSE'")
    if (!file.exists(file) && append)
        append <- FALSE
    Lmismatch  <- ( ! hasLetterAt(as(matches, "DNAStringSet"), gRNA, 
        seq(nchar(gRNA))))
    Lmismatch[!Lmismatch | is.na(Lmismatch)] <- 0
    if (length(dim(Lmismatch)) == 0)
    {
        Lmismatch <- as.data.frame(t(Lmismatch))
    }
    #print(max.mismatch)
    Lmismatch <- Lmismatch[, 1:gRNA.size]
    n.mismatch <- apply(Lmismatch, 1, sum)
    colnames(Lmismatch) <- paste("IsMismatch.pos", 1:gRNA.size, sep = "")
	OffTargetSequence <- toString(as(matches, "DNAStringSet"))
    hits <- data.frame(strand = rep.int(strand, length(matches)),
        chrom = rep.int(seqname, length(matches)), 
        chromStart = start(matches), chromEnd = end(matches), 
        name = names(matches),  
        gRNAPlusPAM = rep(as.character(gRNA), length(matches)),
        OffTargetSequence = unlist(strsplit(OffTargetSequence, ", ")), 
        n.mismatch = n.mismatch,
        chrom.len = rep(chrom.len, length(matches))
    )
    hits <- cbind(Lmismatch,hits)
    hits <- hits[hits$n.mismatch <= max.mismatch,]
    PAM <- gsub("N", "[A|C|G|T]", PAM)
    if (dim(hits)[1] >0)
    {
        containPAM <- unlist(lapply(1:dim(hits)[1], function(i) {
            pos.plus = regexpr(PAM, as.character(hits[i, ]$OffTargetSequence), 
                perl = TRUE)[1]
            if (pos.plus > 0) {
                1
            }
            else { 0 }
        }))
        hits <- hits[containPAM == 1,]
        if (dim(hits)[1] > 0)
        {
            old.start <- hits$chromStart[as.character(hits$strand) == "-"]
            old.end <- hits$chromEnd[as.character(hits$strand) == "-"]
            hits$chromStart[as.character(hits$strand) == "-"] <- (hits$chrom.len 
                - old.end + 1)
            hits$chromEnd[as.character(hits$strand) == "-"] <- (hits$chrom.len 
                - old.start + 1)
            forViewInUCSC <- hits$chrom
            score <- rep(100, dim(hits)[1])
            hits <- hits[, -grep("chrom.len", colnames(hits))]
            hits <- cbind(hits, forViewInUCSC, score)
            hits$forViewInUCSC <- paste(paste(hits$chrom, hits$chromStart, 
                sep = ":"), hits$chromEnd, sep = "-")
            write.table(hits, file = file, append = append, quote = FALSE, 
                sep = "\t", row.names = FALSE, col.names = ! append)
        }
    }
}
