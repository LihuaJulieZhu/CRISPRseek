.mismatches_as_IntegerList <- function(mismatches)
{
    arr_ind <- which(mismatches != 0, arr.ind=TRUE)
    ind_row <- unname(arr_ind[ , "row"])
    ind_col <- unname(arr_ind[ , "col"])
    oo <- orderIntegerPairs(ind_row, ind_col)
    ind_row <- ind_row[oo]
    ind_col <- ind_col[oo]
    partitioning <- PartitioningByEnd(ind_row, NG=nrow(mismatches))
    relist(ind_col, partitioning)
}


#' Build feature vectors
#' 
#' Build feature vectors for calculating scores of off targets
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param hits A Data frame generated from searchHits, which contains
#' \itemize{
#' \item{IsMismatch.posX} - {Indicator variable indicating whether this position X is
#' mismatch or not, (1 means yes and 0 means not). X takes on values from 1 to gRNA.size, representing 
#' all positions in the guide RNA (gRNA).}
#' \item{strand} - {strand of the
#' off target, + for plus and - for minus strand}
#' \item{chrom} - {chromosome of the off
#' target}
#' \item{chromStart} - {start position of the off target}
#' \item{chromEnd} - {end
#' position of the off target}
#' \item{name} - {gRNA name}
#' \item{gRNAPlusPAM} - {gRNA sequence
#' with PAM sequence concatenated}
#' \item{OffTargetSequence} - {the genomic sequence of
#' the off target}
#' \item{n.mismatch} - {number of mismatches between the off target and
#' the gRNA}
#' \item{forViewInUCSC} - {string for viewing in UCSC genome browser, e.g.,
#' chr14:31665685-31665707}
#' \item{score} - {Set to 100, and will be calculated in
#' getOfftargetScore}
#' }
#' @param gRNA.size gRNA size. The default is 20
#' @param canonical.PAM Canonical PAM. The default is NGG for spCas9, TTTN for Cpf1
#' @param subPAM.position The start and end positions of the sub PAM to fetch.
#' Default to 22 and 23 for SP with 20bp gRNA and NGG as preferred PAM
#' @param PAM.size Size of PAM, default to 3 for spCas9, 4 for Cpf1
#' @param PAM.location PAM location relative to gRNA. For example, default to
#' 3prime for spCas9 PAM.  Please set to 5prime for cpf1 PAM since it's PAM is
#' located on the 5 prime end
#' @return A data frame with hits plus features used for calculating scores and
#' for generating report, including 
#' \itemize{
#' \item{IsMismatch.posX} - {Indicator variable
#' indicating whether this position X is mismatch or not, (1 means yes and 0
#' means not, X = 1 - gRNA.size), representing all positions in the gRNA}
#' \item{strand} - {strand of the off target, + for plus and - for minus strand}
#' \item{chrom} - {chromosome of the off target}
#' \item{chromStart} - {start position of the off
#' target}
#' \item{chromEnd} - {end position of the off target}
#' \item{name} - {gRNA name}
#' \item{gRNAPlusPAM} - {gRNA sequence with PAM sequence concatenated}
#' \item{OffTargetSequence} - {the genomic sequence of the off target}
#' \item{n.mismatch} - {number of mismatches between the off target and the gRNA}
#' \item{forViewInUCSC} - {string for viewing in UCSC genome browser, e.g., chr14:31665685-31665707}
#' \item{score} - {score of the off target}
#' \item{mismatche.distance2PAM} - {a comma separated
#' distances of all mismatches to PAM, e.g., 14,11 means one mismatch is 14 bp
#' away from PAM and the other mismatch is 11 bp away from PAM}
#' \item{alignment} - {alignment between gRNA and off target, e.g., ......G..C.......... means
#' that this off target aligns with gRNA except that G and C are mismatches}
#' \item{NGG} - {this off target contains canonical PAM or not, 1 for yes and 0 for no}
#' \item{mean.neighbor.distance.mismatch} - {mean distance between neighboring
#' mismatches}
#' }
#' @note %% ~~further notes~~
#' @author Lihua Julie Zhu
#' @seealso offTargetAnalysis
#' @references %% ~put references to the literature/web site here ~
#' @keywords misc
#' @examples
#' 
#'     hitsFile <-  system.file("extdata", "hits.txt", package = "CRISPRseek")
#'     hits <- read.table(hitsFile, sep= "\t", header = TRUE,
#'         stringsAsFactors = FALSE)
#'     buildFeatureVectorForScoring(hits)
#' @importFrom Biostrings DNAStringSet isMatchingAt substring extractAt complement RNAStringSet DNAString replaceAt
#' @importFrom BiocGenerics sort relist rep.int cbind
#' @importFrom IRanges IRangesList PartitioningByEnd 
#' @importFrom S4Vectors unstrsplit orderIntegerPairs elementNROWS
#' @importFrom methods as
#' @importFrom DelayedArray mean
#' @export 
buildFeatureVectorForScoring <-
    function(hits, gRNA.size = 20, 
             canonical.PAM = "NGG",
             subPAM.position = c(22,23),
             PAM.size = 3, PAM.location = "3prime")
    {
        #hits = read.table(hitsFile, sep = "\t", header=TRUE,
        # stringsAsFactors = FALSE)
        if (dim(hits)[1] == 0)
        {
            stop("Empty hits!")
        }
        subject <- DNAStringSet(as.character(hits$OffTargetSequence))
        if (PAM.location == "3prime")
        {
            isCanonical.PAM <- as.numeric(isMatchingAt(canonical.PAM, subject, 
                                                       at = (gRNA.size + 1), fixed = FALSE))
        }
        else
        {
            isCanonical.PAM <- as.numeric(isMatchingAt(canonical.PAM, subject,
                                                       at = 1, fixed = FALSE))
        }
        PAM <- substring(subject, subPAM.position[1], subPAM.position[2])
        mismatches = hits[, grep("IsMismatch.pos", colnames(hits))]
        mismatch_pos <- .mismatches_as_IntegerList(mismatches)
        #d.nucleotide is the nucleotide that will hybridize to the gRNA
        ### reverse complement of the offtarget sequence
        #r.nucleotide is the gRNA sequence,except T is converted to U)
        
        if (PAM.location == "3prime") 
        {
            at <- IRangesList(start = mismatch_pos, end = mismatch_pos) 
        }
        else
        {
            at.old <- IRangesList(start = mismatch_pos, end = mismatch_pos)
            at <- IRangesList(start = mismatch_pos + PAM.size, 
                              end = mismatch_pos + PAM.size) 
        }
        if (sum(hits$n.mismatch) > 0)
        {
            d.nucleotide <- extractAt(complement(subject), at)
            r.nucleotide <- unlist(extractAt(DNAStringSet(hits$gRNAPlusPAM), at))
            d.nu.r.nu <- paste("r", RNAStringSet(r.nucleotide), ":d", 
                               unlist(d.nucleotide), sep="")
            #### need to assign d.nu.r.nu to the right hits
            arr_ind <- which(mismatches != 0, arr.ind=TRUE)
            ind_row <- sort(unname(arr_ind[ , "row"]))
            partitioning <- PartitioningByEnd(ind_row, NG=nrow(mismatches))
            d.nu.r.nu.2 <- relist(d.nu.r.nu, partitioning)
            d.nu.r.nu <- unstrsplit(as(d.nu.r.nu.2,
                                       "CharacterList"), sep = ",")
        }
        else
        {
            d.nu.r.nu = ""
        }
        
        if (PAM.location == "3prime")
        { 
            mismatch.distance2PAM <- gRNA.size + 1L - mismatch_pos
        }
        else
        {
            mismatch.distance2PAM <- mismatch_pos
        }
        mismatch.distance2PAM <- unstrsplit(as(mismatch.distance2PAM,
                                               "CharacterList"), sep = ",")
        
        alignment <- rep.int(DNAString("."), gRNA.size)
        alignment <- rep.int(DNAStringSet(alignment), nrow(hits))
        if (PAM.location == "3prime")
        {
            alignment <- as.character(replaceAt(alignment, at, extractAt(subject, at)))
        }
        else
        {
            alignment <- as.character(replaceAt(alignment, at.old, extractAt(subject, at)))
        }
        #print(type(mismatch_pos))
        mean.neighbor.distance.mismatch <- mean(diff(mismatch_pos))
        #print(mean.neighbor.distance.mismatch)
        no_neighbor_idx <- elementNROWS(mismatch_pos) <= 1L
       # print(no_neighbor_idx)
        mean.neighbor.distance.mismatch[no_neighbor_idx] <- gRNA.size
        
        features <- cbind(mismatch.distance2PAM, alignment, isCanonical.PAM,
                          mean.neighbor.distance.mismatch, d.nu.r.nu, PAM)
        colnames(features) <- c("mismatch.distance2PAM", "alignment", "NGG", 
                                "mean.neighbor.distance.mismatch", "mismatch.type", "subPAM")
        cbind(hits, features)
    }
