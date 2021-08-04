### for compare2Sequences matches are from matching gRNA portion only


#' Write the hits of sequence search from a sequence to a file
#' 
#' write the hits of sequence search from a sequence instead of BSgenome to a
#' file, internal function used by searchHits
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param gRNA DNAString object with gRNA sequence with PAM appended
#' immediately after,e.g., ACGTACGTACGTACTGACGTCGG with 20bp gRNA sequence plus
#' 3bp PAM sequence CGG
#' @param seqname sequence name as character
#' @param matches XStringViews object storing matched chromosome locations
#' @param strand strand of the match, + for plus strand and - for minus strand
#' @param file file path where the hits is written to
#' @param gRNA.size gRNA size, default 20
#' @param PAM PAM as regular expression for appending to the gRNA, default NGG
#' for SpCas9, change to TTTN for cpf1.
#' @param PAM.pattern PAM as regular expression for filtering the hits, default
#' N[A|G]G$ for spCas9. For cpf1, ^TTTN since it is a 5 prime PAM sequence.
#' @param max.mismatch maximum mismatch allowed within the gRNA (excluding PAM
#' portion) for filtering the hits, default 4
#' @param chrom.len length of the matched chromosome
#' @param append TRUE if append to existing file, false if start a new file
#' @param PAM.location PAM location relative to gRNA. For example, spCas9 PAM
#' is located on the 3 prime while cpf1 PAM is located on the 5 prime
#' @param PAM.size Size of PAM, default 3
#' @param allowed.mismatch.PAM Maximum number of mismatches allowed in the
#' offtargets comparing to the PAM sequence. Default to 1 for NGG PAM
#' @param seqs DNAString object containing a DNA sequence.
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
#' @return results are saved in the file specified by file
#' @note %% ~~further notes~~
#' @author Lihua Julie Zhu
#' @seealso offTargetAnalysis
#' @references
#' http://bioconductor.org/packages/2.8/bioc/vignettes/BSgenome/inst/doc/
#' GenomeSearching.pdf
#' @keywords misc
#' @examples
#' 
#'  if(interactive())
#'  {
#'     gRNAPlusPAM <- DNAString("ACGTACGTACGTACTGACGTCGG")
#'     x <- DNAString("AAGCGCGATATGACGTACGTACGTACTGACGTCGG")
#'     chrom.len <- nchar(as.character(x))
#'     m <- matchPattern(gRNAPlusPAM, x)
#'     names(m) <- "testing"
#'     writeHits(gRNA = gRNAPlusPAM, seqname = "chr1", 
#'         matches = m, strand = "+", file = "exampleWriteHits.txt", 
#'         chrom.len = chrom.len, append = FALSE, seqs = x)
#'  }
#' @importFrom Biostrings hasLetterAt DNAStringSet neditAt DNAString matchPattern
#' @importFrom BiocGenerics unlist cbind rep.int lapply table start end
#' @importFrom methods as
#' @importFrom seqinr s2c
#' @importFrom utils write.table
#' @export  
writeHits <-
    function (gRNA, seqname, matches, strand, file, gRNA.size = 20L, 
        PAM = "NGG", PAM.pattern = "N[A|G]G$", max.mismatch = 4L, 
        chrom.len, append = FALSE,
        PAM.location = "3prime", PAM.size = 3L,
        allowed.mismatch.PAM = 1L,
        seqs,
        baseEditing = FALSE, targetBase = "C", editingWindow = 4:8) 
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
    if (PAM.location == "3prime")
    {
        Lmismatch <- Lmismatch[, 1:gRNA.size]
    }
    else if (dim(Lmismatch)[2] == (gRNA.size + PAM.size))
    {
        start.pos <- PAM.size + 1
        end.pos <- PAM.size + gRNA.size
        Lmismatch <- Lmismatch[, start.pos:end.pos]
    }
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
######### negative strand uses sequences from reverseComplement mapping already
    if(strand == "-")
    {
        new.start1 <- chrom.len - old.end + 1
        new.end1 <- chrom.len -old.start + 1
        if(PAM.location == "3prime")
            new.start1 <- new.start1 - PAM.size
        else
            new.end1 <- new.end1 + PAM.size
    }
    if (PAM.location == "3prime")
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
####### sequence fetch is the same for plus and minus strand
####### because, revcomplement of the sequences (seqs) are used for minus strand
    OffTargetSequence <- substring(as.character(seqs), starts, ends)
###### coordinate needs to be changed for minus strand
    if (strand == "-")
    {
        starts <- unlist(apply(cbind(new.start1,1), 1, max))
        ends <- unlist(apply(cbind(new.end1, chrom.len), 1,min))
    } 
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
    hits
}
