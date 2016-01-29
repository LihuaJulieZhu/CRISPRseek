.mismatches_as_IntegerList <- function(mismatches)
{
    arr_ind <- which(mismatches != 0, arr.ind=TRUE)
    ind_row <- unname(arr_ind[ , "row"])
    ind_col <- unname(arr_ind[ , "col"])
    oo <- S4Vectors:::orderIntegerPairs(ind_row, ind_col)
    ind_row <- ind_row[oo]
    ind_col <- ind_col[oo]
    partitioning <- PartitioningByEnd(ind_row, NG=nrow(mismatches))
    relist(ind_col, partitioning)
}

buildFeatureVectorForScoring <-
    function(hits, gRNA.size = 20, 
    canonical.PAM = "NGG")
{
    #hits = read.table(hitsFile, sep = "\t", header=TRUE,
    # stringsAsFactors = FALSE)
    if (dim(hits)[1] == 0)
    {
        stop("Empty hits!")
    }
    subject <- DNAStringSet(as.character(hits$OffTargetSequence))
    isCanonical.PAM <- as.numeric(isMatchingAt(canonical.PAM, subject, 
        at = (gRNA.size + 1), fixed = FALSE))

    mismatches = hits[, grep("IsMismatch.pos", colnames(hits))]
    mismatch_pos <- .mismatches_as_IntegerList(mismatches)

    mismatch.distance2PAM <- gRNA.size + 1L - mismatch_pos
    mismatch.distance2PAM <- unstrsplit(as(mismatch.distance2PAM,
                                           "CharacterList"),
                                        sep=",")

    alignment <- rep.int(DNAString("."), gRNA.size)
    alignment <- rep.int(DNAStringSet(alignment), nrow(hits))
    at <- IRangesList(start=mismatch_pos, end=mismatch_pos)
    alignment <- as.character(replaceAt(alignment, at, extractAt(subject, at)))

    mean.neighbor.distance.mismatch <- mean(diff(mismatch_pos))
    no_neighbor_idx <- elementNROWS(mismatch_pos) <= 1L
    mean.neighbor.distance.mismatch[no_neighbor_idx] <- gRNA.size

    features <- cbind(mismatch.distance2PAM, alignment, isCanonical.PAM,
        mean.neighbor.distance.mismatch)
    colnames(features) <- c("mismatch.distance2PAM", "alignment", "NGG", 
        "mean.neighbor.distance.mismatch")
    cbind(hits, features)
}
