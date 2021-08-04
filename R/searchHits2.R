.preprocess_me2 <- function(gRNAs, max.mismatch)
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

.searchHitsInOneSeq2 <- function(gRNAs, seq, seqname, PAM, PAM.pattern, PAM.size,
                                max.mismatch, outfile, allowed.mismatch.PAM,
                                PAM.location = "3prime", BSgenomeName,
                                baseEditing = FALSE, targetBase = "C", editingWindow = 4:8)
{
    if (.preprocess_me2(gRNAs, max.mismatch)) {
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
            writeHits2(gRNA = pattern, seqname = seqname, 
                matches = plus_matches, strand = "+",
                file = outfile, gRNA.size = length(pattern),
                PAM = PAM, PAM.pattern = PAM.pattern, 
                max.mismatch = max.mismatch,
                chrom.len = length(seq), append = TRUE, 
                PAM.location = PAM.location,
                PAM.size = PAM.size,
                allowed.mismatch.PAM = allowed.mismatch.PAM,
                BSgenomeName = BSgenomeName,
                baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow)
        }
        if (reverseComplement(pattern) != pattern) {
            minus_matches <- Views(revseq, all_minus_macthes[[i]])
            if (length(minus_matches) > 0) {
                names(minus_matches) <- rep.int(patternID,
                    length(minus_matches))
                writeHits2(gRNA = pattern, seqname = seqname,
                    matches = minus_matches, strand = "-",
                    file = outfile, gRNA.size = length(pattern),
                    PAM = PAM, PAM.pattern = PAM.pattern, 
                    max.mismatch = max.mismatch,
                    chrom.len = length(seq), append = TRUE, 
                    PAM.location = PAM.location, PAM.size = PAM.size,
                    allowed.mismatch.PAM = allowed.mismatch.PAM,
                    BSgenomeName = BSgenomeName,
                    baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow) 
            }
        }
    }
}



#' Search for off targets
#' 
#' Search for off targets for given gRNAs, BSgenome and maximum mismatches
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param gRNAs DNAStringSet object containing a set of gRNAs. Please note the
#' sequences must contain PAM appended after gRNAs, e.g.,
#' ATCGAAATTCGAGCCAATCCCGG where ATCGAAATTCGAGCCAATCC is the gRNA and CGG is
#' the PAM
#' @param BSgenomeName BSgenome object. Please refer to available.genomes in
#' BSgenome package. For example,
#' \itemize{
#' \item{BSgenome.Hsapiens.UCSC.hg19} - for hg19,
#' \item{BSgenome.Mmusculus.UCSC.mm10} - for mm10
#' \item{BSgenome.Celegans.UCSC.ce6} - for ce6
#' \item{BSgenome.Rnorvegicus.UCSC.rn5} - for rn5
#' \item{BSgenome.Drerio.UCSC.danRer7} - for Zv9
#' \item{BSgenome.Dmelanogaster.UCSC.dm3} - for dm3
#' }
#' @param chromToSearch Specify the chromosome to search, default to all,
#' meaning search all chromosomes. For example, chrX indicates searching for
#' matching in chromosome X only
#' @param chromToExclude Specify the chromosome not to search, default to none,
#' meaning to search chromosomes specified by chromToSearch. For example, to
#' exclude haplotype blocks from offtarget search in hg19, set chromToExclude
#' to c(""chr17_ctg5_hap1","chr4_ctg9_hap1", "chr6_apd_hap1", "chr6_cox_hap2",
#' "chr6_dbb_hap3", "chr6_mann_hap4", "chr6_mcf_hap5","chr6_qbl_hap6",
#' "chr6_ssto_hap7")
#' @param max.mismatch Maximum mismatch allowed in off target search, default
#' 3. Warning: will be considerably slower if it is set to greater than 3
#' @param PAM.size Size of PAM, default 3
#' @param gRNA.size Size of gRNA, default 20
#' @param PAM Regular expression of protospacer-adjacent motif (PAM), default
#' NGG for spCas9. For cpf1, ^TTTN
#' @param PAM.pattern Regular expression of PAM, default N[A|G]G$ for spCas9.
#' For cpf1, ^TTTN since it is a 5 prime PAM sequence
#' @param allowed.mismatch.PAM Number of degenerative bases in the PAM
#' sequence, default to 1 for N[A|G]G PAM
#' @param PAM.location PAM location relative to gRNA. For example, spCas9 PAM
#' is located on the 3 prime while cpf1 PAM is located on the 5 prime
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
#' whether this position X is mismatch or not, (1 means yes and 0 means not). X takes on values from 1 to gRNA.size, 
#' representing all positions in the guide RNA (gRNA). }
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
#'     all.gRNAs <- findgRNAs(inputFilePath = 
#'         system.file("extdata", "inputseq.fa", package = "CRISPRseek"),
#'         pairOutputFile = "pairedgRNAs.xls",
#' 	findPairedgRNAOnly = TRUE)
#' 
#'     library("BSgenome.Hsapiens.UCSC.hg19")
#'     ### for speed reason, use max.mismatch = 0 for finding all targets with 
#'     ### all variants of PAM
#'     hits <- searchHits2(all.gRNAs[1], BSgenomeName = Hsapiens,
#'         max.mismatch = 0, chromToSearch = "chrX")
#'     colnames(hits)
#' 
#'     ### test PAM located at 5 prime
#'     all.gRNAs <- findgRNAs(inputFilePath = 
#'              system.file("extdata", "inputseq.fa", package = "CRISPRseek"),
#'              pairOutputFile = "pairedgRNAs.xls",
#'              findPairedgRNAOnly = FALSE,
#'              PAM = "TGT", PAM.location = "5prime")
#'      
#'     library("BSgenome.Hsapiens.UCSC.hg19")
#'          ### for speed reason, use max.mismatch = 0 for finding all targets with 
#'          ### all variants of PAM
#'     hits <- searchHits2(all.gRNAs[1], BSgenomeName = Hsapiens, PAM.size = 3,
#'         max.mismatch = 0, chromToSearch = "chrX", PAM.location = "5prime",
#'         PAM = "NGG", 
#'         PAM.pattern = "^T[A|G]N", allowed.mismatch.PAM = 2)
#'     colnames(hits)
#' @importFrom IRanges width Views
#' @importFrom Biostrings alphabetFrequency matchPDict reverseComplement
#' @importFrom BiocGenerics rep.int
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges intersect setdiff
#' @importFrom XVector subseq
#' @importFrom utils read.table
#' @export
searchHits2 <-
    function (gRNAs, BSgenomeName, chromToSearch = "all", chromToExclude = "", 
	max.mismatch = 3, PAM.size = 3, gRNA.size = 20, 
        PAM = "NGG", PAM.pattern = "N[A|G]G$",
        allowed.mismatch.PAM = 1, PAM.location = "3prime", 
        baseEditing = FALSE, targetBase = "C", editingWindow = 4:8)  
{
    if (missing(gRNAs) || class(gRNAs) != "DNAStringSet") {
        stop("gRNAs is required as a DNAStringSet object!")
    }
    if (missing(BSgenomeName) || class(BSgenomeName) != "BSgenome") {
        stop("BSgenomeName is required as BSgenome object!")
    }
    outfile <- tempfile(tmpdir = getwd())
    max.mismatch <- max.mismatch
    seqnames <- seqnames(BSgenomeName)
    if (chromToSearch != "all")
        seqnames <- intersect(seqnames, chromToSearch)
    if (length(chromToExclude) >0)
	seqnames <- setdiff(seqnames, chromToExclude)
    append <- FALSE
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
    for (seqname in seqnames) {
        cat(">>> Finding all hits in sequence", seqname, "...\n")
        subject <- BSgenomeName[[seqname]]
        .searchHitsInOneSeq2(gRNAs = gRNAs, seq = subject, 
             seqname = seqname, PAM = PAM, 
             PAM.pattern = PAM.pattern, PAM.size = PAM.size,
             max.mismatch = max.mismatch, outfile, 
             allowed.mismatch.PAM = allowed.mismatch.PAM,
             PAM.location = PAM.location,
             BSgenomeName = BSgenomeName,
             baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow)
        cat(">>> DONE searching\n")
    }
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
