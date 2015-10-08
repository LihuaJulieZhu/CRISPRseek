buildFeatureVectorForScoring <-
    function(hits, gRNA.size = 20, canonical.PAM = "NGG")
{
    #hits = read.table(hitsFile, sep = "\t", header=TRUE,
    # stringsAsFactors = FALSE)
    if (dim(hits)[1] == 0)
    {
        stop("Empty hits!")
    }
    subject <- DNAStringSet(as.character(hits$OffTargetSequence))
    pattern <- DNAString(as.character(hits$gRNAPlusPAM[1]))
    isCanonical.PAM <- as.numeric(isMatchingAt(canonical.PAM, subject, 
        at = (gRNA.size + 1), fixed = FALSE))
    #type.mismatch = matrix(nrow=dim(hits)[1], ncol=gRNA.size)
    mismatch.pos = hits[, grep("IsMismatch.pos", colnames(hits))]
    mismatch.distance2PAM <- apply(mismatch.pos, 1, function(i) { 
        paste(gRNA.size + 1 - which(i == 1), collapse = ",")
    })
    n.cores <- detectCores() - 1
    cl <- makeCluster(n.cores)
    clusterExport(cl, varlist = c("mismatch.pos", "gRNA.size",
        "rev", "subject", "hits",  "subseq", "mismatch.distance2PAM",
        "strsplit"),
        envir = environment())
    alignment <- unlist(parLapply(cl, 1:dim(mismatch.pos)[1], function(i) {
        temp <- rep(".", gRNA.size)
        ind <- which(mismatch.pos[i,] == 1)
        for (j in ind)
            temp[j] = as.character(subseq(subject[i], start = j, width = 1))
        paste(temp, collapse = "")
    }))
    mean.neighbor.distance.mismatch <- unlist(parLapply(cl, 1:dim(hits)[1],
       function(i) 
    {
        positions <- rev(as.numeric(unlist(strsplit(mismatch.distance2PAM[i],
            split = ","))))
        n.mismatch <- hits$n.mismatch[i]
        if (n.mismatch > 1)
            mean(positions[2:n.mismatch] - positions[1:(n.mismatch-1)])
        else
            gRNA.size
    }))
    stopCluster(cl)
    features <- cbind(mismatch.distance2PAM, alignment, isCanonical.PAM,
        mean.neighbor.distance.mismatch)
    colnames(features) <- c("mismatch.distance2PAM", "alignment", "NGG", 
        "mean.neighbor.distance.mismatch")
    cbind(hits, features)
}
