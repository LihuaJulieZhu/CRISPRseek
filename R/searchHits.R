####### for compare2Sequences
.searchHitsInOneSeq3 <- function(gRNAs, seqs, seqname, 
    PAM = "NGG", PAM.pattern = "NNN$", PAM.size = 3,
    max.mismatch = 2, outfile, allowed.mismatch.PAM = 2,
    PAM.location = "3prime", gRNA.size = 20,
    baseEditing = FALSE, targetBase = "C", editingWindow = 4:8)
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
                seqs = seqs,
		baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow)
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
                    seqs = revseq, 
                    baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow)
            }
        }
    }
}



#' Search for off targets in a sequence as DNAString
#' 
#' Search for off targets for given gRNAs, sequence and maximum mismatches
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param gRNAs DNAStringSet object containing a set of gRNAs. Please note the
#' sequences must contain PAM appended after gRNAs, e.g.,
#' ATCGAAATTCGAGCCAATCCCGG where ATCGAAATTCGAGCCAATCC is the gRNA and CGG is
#' the PAM
#' @param seqs DNAString object containing a DNA sequence.
#' @param seqname Specify the name of the sequence
#' @param max.mismatch Maximum mismatch allowed in off target search, default
#' 3. Warning: will be considerably slower if it is set to greater than 3
#' @param PAM.size Size of PAM, default 3
#' @param gRNA.size Size of gRNA, default 20
#' @param PAM PAM as regular expression for appending to the gRNA, default NGG
#' for SpCas9, change to TTTN for cpf1.
#' @param PAM.pattern Regular expression of PAM, default N[A|G]G$ for spCas9.
#' For cpf1, ^TTTN since it is a 5 prime PAM sequence
#' @param allowed.mismatch.PAM Maximum number of mismatches allowed in the
#' offtargets comparing to the PAM sequence. Default to 2 for NGG PAM
#' @param PAM.location PAM location relative to gRNA. For example, spCas9 PAM
#' is located on the 3 prime while cpf1 PAM is located on the 5 prime
#' @param outfile File path to temporarily store the search results
#' @param baseEditing Indicate whether to design gRNAs for base editing.
#' Default to FALSE If TRUE, please set baseEditing = TRUE, targetBase and
#' editingWidow accordingly.
#' @param targetBase Applicable only when baseEditing is set to TRUE. It is
#' used to indicate the target base for base editing systems, default to C for
#' converting C to T in the CBE system. Please change it to A if you intend to
#' use the ABE system.
#' @param editingWindow Applicable only when baseEditing is set to TRUE. It is
#' used to indicate the effective editing window to consider for the offtargets
#' search only, default to 4 to 8 which is for the original CBE system. Please
#' change it accordingly if the system you use have a different editing window,
#' or you would like to include offtargets with the target base in a larger
#' editing window.
#' @return a data frame contains 
#' \itemize{
#' \item{IsMismatch.posX} - {indicator variable indicating
#' whether this position X is mismatch or not, (1 means yes and 0 means not). X takes on values from 1 to gRNA.size, representing all positions in the guide RNA (gRNA). }
#' \item{strand} - {strand of
#' the match, + for plus and - for minus strand}
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
#' \item{score} - {set to 100, and will be updated in
#' getOfftargetScore}
#' }
#' @note %% ~~further notes~~
#' @author Lihua Julie Zhu
#' @seealso offTargetAnalysis
#' @references %% ~put references to the literature/web site here ~
#' @keywords misc
#' @examples
#' 
#'  all.gRNAs <- findgRNAs(inputFilePath =
#'        system.file("extdata", "inputseq.fa", package = "CRISPRseek"),
#'        pairOutputFile = "pairedgRNAs.xls",
#'        findPairedgRNAOnly = TRUE)

#'   hits <- searchHits(all.gRNAs[1], 
#'       seqs = DNAString(
#'            "TAATATTTTAAAATCGGTGACGTGGGCCCAAAACGAGTGCAGTTCCAAAGGCACCCACCTGTGGCAG"), 
#'       seqname = "myseq", max.mismatch = 10, outfile = "test_searchHits")

#'   colnames(hits)

    ### test PAM located at 5 prime
#'   all.gRNAs <- findgRNAs(inputFilePath =
#'           DNAStringSet(
#'               "TAATATTTTAAAATCGGTGACGTGGGCCCAAAACGAGTGCAGTTCCAAAGGCACCCACCTGTGGCAG"),
#'           pairOutputFile = "pairedgRNAs.xls",
#'           findPairedgRNAOnly = FALSE,
#'           PAM = "TTTN", PAM.location = "5prime")

#'    hits <- searchHits(all.gRNAs[1], seqs = DNAString(
#'        "TAATATTTTAAAATCGGTGACGTGGGCCCAAAACGAGTGCAGTTCCAAAGGCACCCACCTGTGGCAG"), 
#'        seqname = "myseq", 
#'        max.mismatch = 0, 
#'        outfile = "test_searchHits", PAM.location = "5prime",
#'        PAM.pattern = "^T[A|T]NN", allowed.mismatch.PAM = 0, PAM = "TTTN")
#'    colnames(hits)

#' @importFrom Biostrings PDict matchPDict reverseComplement
#' @importFrom IRanges Views width
#' @importFrom BiocGenerics rep.int
#' @importFrom utils read.table
#' @importFrom XVector subseq
#' @export
searchHits <-
    function (gRNAs,seqs,seqname,
	max.mismatch = 3, PAM.size = 3, gRNA.size = 20, 
        PAM = "NGG", PAM.pattern = "NNN$",
        allowed.mismatch.PAM = 2, PAM.location = "3prime",
        outfile,
        baseEditing = FALSE, targetBase = "C", editingWindow = 4:8) 
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
             PAM.location = PAM.location,
             baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow)
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
