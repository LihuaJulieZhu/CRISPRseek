.mismatches_as_IntegerList <- function(mismatches)
{
    arr_ind <- which(mismatches != 0, arr.ind=TRUE)
    ind_row <- unname(arr_ind[ , "row"])
    ind_col <- unname(arr_ind[ , "col"])
    oo <- S4Vectors::orderIntegerPairs(ind_row, ind_col)
    ind_row <- ind_row[oo]
    ind_col <- ind_col[oo]
    partitioning <- PartitioningByEnd(ind_row, NG=nrow(mismatches))
    relist(ind_col, partitioning)
}

#' @importFrom Biostrings DNAStringSet isMatchingAt substring extractAt complement RNAStringSet DNAString replaceAt
#' @importFrom BiocGenerics grep sort relist rep.int cbind
#' @importFrom IRanges IRangesList PartitioningByEnd
#' @importFrom S4Vectors unstrsplit orderIntegerPairs elementNROWS
#' @importFrom methods as
#' @importFrom DelayedArray mean
#'
#'
#' @example
#' alns <- readRDS("~/DropboxUmass/Bioconductor/Trunk/GUIDEseq/inst/extdata/alns.RDS")
#' temp <- buildFeatureVectorForScoringBulge(alns)
#' expect_equal(temp$mismatch.type[1], "rU:dC,rU:dG,rA:dA")
#' expect_equal(as.numeric(temp$mean.neighbor.distance.mismatch), c(9,20))
#' expect_equal(temp$gRNA.deletion[1], "G")
#'
buildFeatureVectorForScoringBulge <-
    function(alns, gRNA.size = 20,
    canonical.PAM = "NGG", insertion.symbol = "^",
     PAM.size = 3, PAM.location = "3prime")
{
    #hits = read.table(hitsFile, sep = "\t", header=TRUE,
    # stringsAsFactors = FALSE)
    if (PAM.location != "3prime")
        stop("Bulge scoring for 5prime PAM not implemented yet")
    if (dim(alns)[1] == 0)
    {
        stop("Empty hits!")
    }
    #subject <- DNAStringSet(as.character(alns$offTarget_sequence))

    PAM <- DNAStringSet(unlist(alns$PAM.sequence))
    isCanonical.PAM <- as.numeric(isMatchingAt(canonical.PAM,
                                               PAM,
                                               at = 1, fixed = FALSE))
    mismatch_pos <- alns$pos.mismatch
    insertion_pos <- alns$pos.insertion
    deletion_pos <- alns$pos.deletion

    mismatches <- matrix(rep(0, gRNA.size * length(PAM)),
                         ncol = gRNA.size, nrow = length(PAM))
    colnames(mismatches) <- paste0("IsMismatch.pos", 1:gRNA.size)
    for (i in 1:nrow(mismatches))
        mismatches[i, mismatch_pos[[i]]] <- 1

    insertions <- matrix(rep(0, gRNA.size * length(PAM)),
                         ncol = gRNA.size, nrow = length(PAM))
    colnames(insertions) <- paste0("IsInsertion.pos", 1:gRNA.size)
    for (i in 1:nrow(insertions))
      insertions[i, insertion_pos[[i]]] <- 1

    deletions <- matrix(rep(0, gRNA.size * length(PAM)),
                         ncol = gRNA.size, nrow = length(PAM))
    colnames(deletions) <- paste0("IsDeletion.pos", 1:gRNA.size)
    for (i in 1:nrow(deletions))
      deletions[i, deletion_pos[[i]]] <- 1

    #d.nucleotide is the nucleotide that will hybridize to the gRNA
    ### reverse complement of the offtarget sequence
    #r.nucleotide is the gRNA sequence,except T is converted to U)

    at <- IRangesList(start=mismatch_pos, end=mismatch_pos)
    at.ins <- IRangesList(start=insertion_pos, end=insertion_pos)
    at.del <- IRangesList(start=deletion_pos, end=deletion_pos)

    if (insertion.symbol == "lowerCase")
    {
        insertion.symbol =="[acgt]"
        OT <- gsub(".", "A", gsub(insertion.symbol, "",
                   alns$guideAlignment2OffTarget,
                   fixed = FALSE), fixed = TRUE)
    }
    else
    {
        OT <- gsub(".", "A", gsub(insertion.symbol, "",
                   alns$guideAlignment2OffTarget,
                   fixed = TRUE), fixed = TRUE)
    }

    r.ins <- as.character(unlist(extractAt(
      DNAStringSet(unlist(alns$gRNAPlusPAM)), at.ins)))
    r.del <- as.character(unlist(extractAt(
      DNAStringSet(unlist(alns$offTarget_sequence)), at.del)))

    d.nucleotide <- extractAt(complement(DNAStringSet(OT)), at)
    r.nucleotide <- as.character(unlist(extractAt(
        DNAStringSet(alns$gRNAPlusPAM), at)))
    r.nucleotide[r.nucleotide == "T"] <- "U"
    d.nu.r.nu <- paste("r", r.nucleotide, ":d",
      unlist(d.nucleotide), sep="")
    #### need to assign d.nu.r.nu to the right alignment
    arr_ind <- which(mismatches != 0, arr.ind=TRUE)
    ind_row <- sort(unname(arr_ind[ , "row"]))
    partitioning <- PartitioningByEnd(ind_row, NG=nrow(mismatches))
    d.nu.r.nu.2 <- relist(d.nu.r.nu, partitioning)
    d.nu.r.nu <- unstrsplit(as(d.nu.r.nu.2,
                   "CharacterList"), sep = ",")

    #### need to assign insertion and deletion to the right alignment
    arr_ind <- which(insertions != 0, arr.ind=TRUE)
    ind_row <- sort(unname(arr_ind[ , "row"]))
    partitioning <- PartitioningByEnd(ind_row, NG=nrow(insertions))
    r.ins.2 <- relist(r.ins, partitioning)
    r.ins <- unstrsplit(as(r.ins.2,
                               "CharacterList"), sep = ",")

    arr_ind <- which(deletions != 0, arr.ind=TRUE)
    ind_row <- sort(unname(arr_ind[ , "row"]))
    partitioning <- PartitioningByEnd(ind_row, NG=nrow(deletions))
    r.del.2 <- relist(r.del, partitioning)
    r.del <- unstrsplit(as(r.del.2,
                           "CharacterList"), sep = ",")

    mismatch.distance2PAM <-  unlist(alns$mismatch.distance2PAM)

    # alignment <- rep.int(DNAString("."), gRNA.size)
    # alignment <- rep.int(DNAStringSet(alignment), nrow(hits))
    # alignment <- as.character(replaceAt(alignment, at, extractAt(subject, at)))

    mean.neighbor.distance.mismatch <- unlist(lapply(1:nrow(mismatches),
                                              function(i) {
                                                  mean(diff(mismatch_pos[[i]]))
                                              }))
    no_neighbor_idx <- elementNROWS(mismatch_pos) <= 1L
    mean.neighbor.distance.mismatch[no_neighbor_idx] <- gRNA.size

    subPAM <- substr(alns$PAM.sequence, 2,3)
    features <- cbind(isCanonical.PAM,
        mean.neighbor.distance.mismatch, d.nu.r.nu, subPAM, r.ins, r.del)
    colnames(alns)[colnames(alns) == "guideAlignment2OffTarget"] <- "alignment"

    colnames(features) <- c("NGG",
        "mean.neighbor.distance.mismatch", "mismatch.type", "subPAM",
        "gRNA.insertion",
        "gRNA.deletion")
    cbind(alns, mismatches, insertions, deletions, features)
}
