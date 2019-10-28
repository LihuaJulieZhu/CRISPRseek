writeHits2 <-
    function (gRNA, seqname, matches, strand, file, gRNA.size = 20L, 
        PAM = "NGG", PAM.pattern = "N[A|G]G$", max.mismatch = 4L, 
        chrom.len, append = FALSE,
        PAM.location = "3prime", PAM.size = 3L,
        allowed.mismatch.PAM = 1L,
        BSgenomeName, baseEditing = FALSE, targetBase = "C", editingWindow = 4:8) 
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
    #if (PAM.location == "3prime")
    #{
    #    Lmismatch <- Lmismatch[, 1:gRNA.size]
    #}
    #else if (dim(Lmismatch)[2] == (gRNA.size + PAM.size))
    #{
    #    start.pos <- PAM.size + 1
    #    end.pos <- PAM.size + gRNA.size
    #    Lmismatch <- Lmismatch[, start.pos:end.pos]
    #}
    n.mismatch <- apply(Lmismatch, 1, sum)
    colnames(Lmismatch) <- paste("IsMismatch.pos", 1:gRNA.size, sep = "")
    
#############################  
### for changing the definition 
### of allowed.mismatch.PAM
### to mismatch to canonical PAM
#############################
    if(PAM.location == "3prime")
        gRNAplusPAM <- paste(as.character(gRNA), 
            as.character(PAM), sep="")
    else
        gRNAplusPAM <- paste(as.character(PAM),
            as.character(gRNA), sep="")
    old.start <- start(matches)
    old.end <- end(matches)
    if(strand == "-")
    {
        new.start <- chrom.len - old.end + 1
        new.end <- chrom.len -old.start + 1
        if(PAM.location == "3prime")
            new.start <- new.start - PAM.size
        else
            new.end <- new.end + PAM.size
    }
    else if (PAM.location == "3prime")
    {
        new.start <- old.start
        new.end <- old.end + PAM.size
    }
    else
    {
        new.start <- old.start - PAM.size
        new.end <- old.end
    }
    starts <- unlist(apply(cbind(new.start,1), 1, max))
    ends <- unlist(apply(cbind(new.end, chrom.len), 1,min))
    OffTargetSequence <- unlist(strsplit(toString(BSgenome::getSeq(BSgenomeName,
            names = seqname, start = starts,
            end = ends, strand = strand, width = NA, as.character = FALSE)), ", "))
#############################
    #OffTargetSequence2 <- unlist(strsplit(toString(as(matches, "DNAStringSet")), ", "))
    hits <- data.frame(strand = rep.int(strand, length(matches)),
        chrom = rep.int(seqname, length(matches)), 
        chromStart = starts, chromEnd = ends, 
        name = names(matches),  
        gRNAPlusPAM = rep(gRNAplusPAM, length(matches)),
        OffTargetSequence = OffTargetSequence, 
        n.mismatch = n.mismatch,
        chrom.len = rep(chrom.len, length(matches))
    )
    hits <- cbind(Lmismatch,hits)
    hits <- subset(hits, n.mismatch <= max.mismatch & 
         (ends - starts + 1) == (gRNA.size + PAM.size))
 
    PAM.pattern <- translatePattern(PAM.pattern)
    if (dim(hits)[1] >0)
    {
         containPAM <- unlist(lapply(1:dim(hits)[1], function(i) {
             pos.plus = regexpr(PAM.pattern, 
                 as.character(hits[i, ]$OffTargetSequence), perl = TRUE)[1]
             if (pos.plus > 0) {
                 1
             }
             else { 0 }
        }))
        hits <- hits[containPAM == 1,]
        if (dim(hits)[1] > 0)
        {
            if (baseEditing)
            {
          
                n.targetBase <- unlist(lapply(1:dim(hits)[1], function(i) {
                   table(factor(s2c(substring(as.character(hits[i, ]$OffTargetSequence), 
                   min(editingWindow),max(editingWindow))), levels=c(targetBase)))
                }))
                hits <- hits[n.targetBase > 0, ]
            }   
            if (dim(hits)[1] > 0)
            {
              if (PAM.location == "3prime")
              {
                PAM.sequence <- substr(hits$OffTargetSequence,
                    gRNA.size + 1, gRNA.size + PAM.size)
              }
              else
              {
                PAM.sequence <- substr(hits$OffTargetSequence,
                    1,  PAM.size)
              }

              n.PAM.mismatch <- unlist(lapply(DNAStringSet(PAM.sequence), function(i) {
                  neditAt(i, DNAString(PAM), fixed=FALSE)
                  }))
             
              hits <- hits[n.PAM.mismatch <= allowed.mismatch.PAM,]
              if (dim(hits)[1] > 0)
              {
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
    }
}
