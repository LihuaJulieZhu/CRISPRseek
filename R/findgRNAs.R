.getgRNA.cut.sites <- function(subject, subjectname, PAM ="NGG", 
    gRNA.pattern = "", gRNA.size = 20, cut.site = 17,
    PAM.size = 3, calculategRNAEfficacy = TRUE, baseBeforegRNA = 4,
    baseAfterPAM = 3,
    reverse.subject = FALSE,
    PAM.location = "3prime",
    baseEditing = FALSE, targetBase = "C", editingWindow = 4:8
    )
{
     seq.len <- nchar(as.character(subject))
     pos.PAMs <- unlist(gregexpr(PAM, subject, perl = TRUE,
         ignore.case = TRUE, fixed = FALSE))
     if (PAM.location == "3prime")
         pos.PAMs <- pos.PAMs[pos.PAMs != -1 & pos.PAMs > gRNA.size]
     else
         pos.PAMs <- pos.PAMs[pos.PAMs != -1]
     if (length(pos.PAMs) > 0)
     {
        if (PAM.location == "3prime")
            starts.gRNA <- pos.PAMs - gRNA.size
        else
            starts.gRNA <- pos.PAMs + PAM.size 
        starts.gRNA <- subset(starts.gRNA, starts.gRNA <= (seq.len - gRNA.size + 1)) 
        starts.gRNA <- subset(starts.gRNA, starts.gRNA > 0)
        ends.gRNA <- starts.gRNA + gRNA.size - 1
        if (gRNA.pattern != "" && length(starts.gRNA))
        {
            gRNA.seqs <- as.character(Views(subject, 
               start = starts.gRNA,
               end = ends.gRNA))
            pos.set2 <- unlist(gregexpr(gRNA.pattern, gRNA.seqs,
               perl = TRUE, ignore.case = TRUE, fixed = FALSE))
            pos.PAMs <- pos.PAMs[pos.set2 == 1]
       }
       if (baseEditing && length(pos.PAMs) > 0)
       {
            if (PAM.location == "3prime")
                starts.gRNA <- pos.PAMs - gRNA.size
            else
                starts.gRNA <- pos.PAMs + PAM.size 
            ends.gRNA <- starts.gRNA + gRNA.size - 1
          
            pos.PAMs <- pos.PAMs[ends.gRNA <= length(subject)]
	    starts.gRNA <- starts.gRNA[ends.gRNA <= length(subject)]
            ends.gRNA <- ends.gRNA[ends.gRNA <= length(subject)]
            gRNA.seqs <- as.character(Views(subject,
                 start = starts.gRNA,
                 end = ends.gRNA))
             n.targetBase <- unlist(lapply(1:length(gRNA.seqs), function(i) {
                table(factor(s2c(substring(as.character(gRNA.seqs[i]), min(editingWindow),max(editingWindow))), levels=c(targetBase)))
             }))
            pos.PAMs <- pos.PAMs[n.targetBase == 1]
            starts.gRNA <- starts.gRNA[n.targetBase == 1]
            ends.gRNA <- ends.gRNA[n.targetBase == 1]
       }       
      if (length(pos.PAMs) > 0 && length(starts.gRNA))
      {
         if (PAM.location == "3prime")
            seq <- as.character(Views(subject,
                start = starts.gRNA,
                end = ends.gRNA + PAM.size))
         else
            seq <- as.character(Views(subject,
                start = starts.gRNA - PAM.size,
                end = ends.gRNA))
      
          extendedSequence <- seq
          if (calculategRNAEfficacy)
          {
            extended.starts <- starts.gRNA - baseBeforegRNA
            extended.starts[extended.starts < 1] <- 1
            if (PAM.location == "3prime")
                extended.ends <- ends.gRNA + PAM.size + baseAfterPAM
            else
                extended.ends <- starts.gRNA + baseAfterPAM - 1
            extended.ends[extended.ends > length(subject)] <- length(subject)
            extendedSequence <- as.character(Views(subject, start = extended.starts,
                end = extended.ends))
         }
         if (reverse.subject)
           {
              gRNAs.cut <- cbind(seq, paste( subjectname,"_gR",
                (seq.len - (starts.gRNA + cut.site -1 ) + 1), "r", sep = ""),
                (seq.len - starts.gRNA + 1), "-", extendedSequence)
           }
         else
           {
              gRNAs.cut <-
                cbind(seq, paste(subjectname,"_gR", (starts.gRNA + cut.site - 1),
                "f", sep = ""), starts.gRNA, "+", extendedSequence)
           }
      }
      else
      {
 	gRNAs.cut <- "" 
      }
     }
     else
     {
        gRNAs.cut <- ""
     }
        gRNAs.cut
}

## Very inefficient (perform in quadratic time).
.compute_pair_index <- function(plus_start, minus_start, min.gap, max.gap)
{
    plus.index <- minus.index <- integer(0)
    for (j in seq_along(plus_start))
    {
        for (k in seq_along(minus_start))
        {
            d <- plus_start[j] - minus_start[k]
            if (d > min.gap && d <= max.gap)
            {
                plus.index <- c(plus.index, j)
                minus.index <- c(minus.index, k)
            } ### if paired
        } #for k
    }## for j
    list(plus.index, minus.index)
}

## A more efficient version of .compute_pair_index() that performs almost in
## linear time.
.compute_pair_index2 <- function(plus_start, minus_start, min.gap, max.gap)
{
    plus.index <- minus.index <- integer(0)
    n1 <- length(plus_start)
    n2 <- length(minus_start)
    if (n1 == 0L || n2 == 0L)
        return(list(plus.index, minus.index))
    plus_oo <- order(plus_start)
    minus_oo <- order(minus_start)
    i2 <- 1L
    for (i1 in seq_len(n1)) {
        j <- plus_oo[[i1]]
        current_plus_start <- plus_start[[j]]
        min_minus_start <- current_plus_start - max.gap
        max_minus_start <- current_plus_start - min.gap - 1
        while (minus_start[[minus_oo[[i2]]]] < min_minus_start) {
            i2 <- i2 + 1L
            if (i2 > n2)
                break
        }
        if (i2 > n2)
            break
        i <- i2
        while (minus_start[[k <- minus_oo[[i]]]] <= max_minus_start) {
            plus.index <- c(plus.index, j)
            minus.index <- c(minus.index, k)
            i <- i + 1L
            if (i > n2)
                break
        }
    }
    oo <- orderIntegerPairs(plus.index, minus.index)
    list(plus.index[oo], minus.index[oo])
}

## An even more efficient version of .compute_pair_index() based on
## findOverlaps(). Hundreds times faster than .compute_pair_index2().
.compute_pair_index3 <- function(plus_start, minus_start, min.gap, max.gap)
{
    #shift <- floor((min.gap + max.gap) / 2)
    #subject <- IRanges(minus_start + shift, width=1L)
    #maxgap <- ceiling((max.gap - min.gap) / 2)
    subject <- IRanges(minus_start, width = 1L)
    query <- IRanges(plus_start, width=1L)
    hits <- findOverlaps(query, subject, maxgap=max.gap)
    d <- plus_start[queryHits(hits)] - minus_start[subjectHits(hits)]
    hits <- sort(hits[min.gap <= d & d <= max.gap])
    gap <- plus_start[queryHits(hits)] - minus_start[subjectHits(hits)]
    list(queryHits(hits), subjectHits(hits), gap)
}



#' Find potential gRNAs
#' 
#' Find potential gRNAs for an input file containing sequences in fasta format
#' 
#' If users already has a fasta file that contains a set of potential gRNAs,
#' then users can call filergRNAs directly although the easiest way is to call
#' the one-stop-shopping function OffTargetAnalysis with findgRNAs set to
#' FALSE.
#' 
#' @param inputFilePath Sequence input file path or a DNAStringSet object that
#' contains sequences to be searched for potential gRNAs
#' @param baseEditing Indicate whether to design gRNAs for base editing.
#' Default to FALSE If TRUE, please set baseEditing = TRUE, targetBase and
#' editingWidow accordingly.
#' @param targetBase Applicable only when baseEditing is set to TRUE. It is
#' used to indicate the target base for base editing systems, default to C for
#' converting C to T in the CBE system. Please change it to A if you intend to
#' use the ABE system.
#' @param editingWindow Applicable only when baseEditing is set to TRUE. It is
#' used to indicate the effective editing window, default to 4 to 8 which is
#' for the original CBE system. Please change it accordingly if the system you
#' use have a different editing window.
#' @param format Format of the input file, fasta and fastq are supported,
#' default fasta
#' @param PAM protospacer-adjacent motif (PAM) sequence after the gRNA, default
#' NGG
#' @param PAM.size PAM length, default 3
#' @param findPairedgRNAOnly Choose whether to only search for paired gRNAs in
#' such an orientation that the first one is on minus strand called reverse
#' gRNA and the second one is on plus strand called forward gRNA. TRUE or
#' FALSE, default FALSE
#' @param annotatePaired Indicate whether to output paired information, default
#' TRUE
#' @param paired.orientation PAMin orientation means the two adjacent PAMs on
#' the sense and antisense strands face inwards towards each other like N21GG
#' and CCN21 whereas PAMout orientation means they face away from each other
#' like CCN21 and N21GG
#' @param enable.multicore Indicate whether enable parallel processing, default
#' FALSE. For super long sequences with lots of gRNAs, suggest set it to TRUE
#' @param n.cores.max Indicating maximum number of cores to use in multi core
#' mode, i.e., parallel processing, default 6. Please set it to 1 to disable
#' multicore processing for small dataset.
#' @param gRNA.pattern Regular expression or IUPAC Extended Genetic Alphabet to
#' represent gRNA pattern, default is no restriction. To specify that the gRNA
#' must start with GG for example, then set it to ^GG. Please see
#' help(translatePattern) for a list of IUPAC Extended Genetic Alphabet.
#' @param gRNA.size The size of the gRNA, default 20
#' @param overlap.gRNA.positions The required overlap positions of gRNA and
#' restriction enzyme cut site, default 17 and 18. For Cpf1,you may set it to 19 and 23.
#' @param primeEditing Indicate whether to design gRNAs for prime editing.
#' Default to FALSE.  If true, please set PBS.length, RT.template.length,
#' RT.template.pattern, targeted.seq.length.change, bp.after.target.end,
#' target.start, and target.end accordingly
#' @param PBS.length Applicable only when primeEditing is set to TRUE. It is
#' used to specify the number of bases to ouput for primer binding site.
#' @param RT.template.length Applicable only when primeEditing is set to TRUE.
#' It is used to specify the number of bases required for RT template, default
#' to 8 to 18. Please increase the length if the edit is large insertion.  Only
#' gRNAs with calculated RT.template.length falling into the specified range
#' will be in the output. It is calculated as the following. RT.template.length
#' = target.start – cut.start + (target.end - target.start) +
#' targeted.seq.length.change + bp.after.target.end
#' @param RT.template.pattern Applicable only when primeEditing is set to TRUE.
#' It is used to specify the RT template sequence pattern, default to not
#' ending with C according to https://doi.org/10.1038/s41586-019-1711-4
#' @param corrected.seq Applicable only when primeEditing is set to TRUE. It is
#' used to specify the mutated or inserted sequences after successful editing.
#' @param targeted.seq.length.change Applicable only when primeEditing is set
#' to TRUE. It is used to specify the number of targeted sequence length
#' change. Please set it to 0 for base changes, positive numbers for insersion,
#' and negative number for deletion. For example, 10 means that the corrected
#' sequence will have 10bp insertion, -10 means that the corrected sequence
#' will have 10bp deletion, and 0 means only bases have been changed and the
#' sequence length remains the same
#' @param bp.after.target.end Applicable only when primeEditing is set to TRUE.
#' It is used to specify the number of bases to add after the target change end
#' site as part of RT template. Please refer to RT.template.length for how this
#' parameter influences the RT.template.length calculation which is used as a
#' filtering criteria in pregRNA selection.
#' @param target.start Applicable only when primeEditing is set to TRUE. It is
#' used to specify the start location in the input sequence to make changes,
#' which will be used to obtain the RT template sequence. Please also refer to
#' RT.template.length for how this parameter influences the RT.template.length
#' calculation which is used as a filtering criteria in pregRNA selection.
#' @param target.end Applicable only when primeEditing is set to TRUE. It is
#' used to specify the end location in the input sequnence to make changes,
#' which will be used to obtain the RT template sequence. Please also refer to
#' RT.template.length for how this parameter influences the RT.template.length
#' calculation which is used as a filtering criteria in pregRNA selection.
#' @param primeEditingPaired.output Applicable only when primeEditing is set to
#' TRUE. It is used to specify the file path to save pegRNA and the second gRNA
#' with PBS, RT.template, gRNA sequences, default pairedgRNAsForPE.xls
#' @param min.gap Minimum distance between two oppositely oriented gRNAs to be
#' valid paired gRNAs. Default 0
#' @param max.gap Maximum distance between two oppositely oriented gRNAs to be
#' valid paired gRNAs. Default 20
#' @param pairOutputFile The output file for writing paired gRNA information to
#' @param name.prefix The prefix used when assign name to found gRNAs, default
#' gRNA, short for guided RNA.
#' @param baseBeforegRNA Number of bases before gRNA used for calculating gRNA
#' efficiency, default 4 for spCas9 Please note, for PAM located on the 5
#' prime, need to specify the number of bases before the PAM sequence plus PAM
#' size.
#' @param baseAfterPAM Number of bases after PAM used for calculating gRNA
#' efficiency, default 3 for spCas9 Please note, for PAM located on the 5
#' prime, need to include the length of the gRNA plus the extended sequence on
#' the 3 prime
#' @param featureWeightMatrixFile Feature weight matrix file used for
#' calculating gRNA efficiency. By default DoenchNBT2014 weight matrix is used.
#' To use alternative weight matrix file, please input a csv file with first
#' column containing significant features and the second column containing the
#' corresponding weights for the features. Please see Doench et al., 2014 for
#' details.
#' @param calculategRNAEfficacy Default to FALSE, not to calculate gRNA
#' efficacy
#' @param efficacyFile File path to write gRNA efficacies
#' @param PAM.location PAM location relative to gRNA. For example, spCas9 PAM
#' is located on the 3 prime while cpf1 PAM is located on the 5 prime
#' @param rule.set Specify a rule set scoring system for calculating gRNA
#' efficacy. Please note that if specifying DeepCpf1, please specify other
#' parameters accordingly for CRISPR-Cpf1 gRNAs.
#' @param chrom_acc Optional binary variable indicating chromatin accessibility 
#' information with 1 indicating accessible and 0 not accessible.
#' @return DNAStringSet consists of potential gRNAs that can be input to
#' filtergRNAs function directly
#' @note If the input sequence file contains multiple >300 bp sequences,
#' suggest create one input file for each sequence and run the
#' OffTargetAnalysis separately.
#' @author Lihua Julie Zhu
#' @seealso offTargetAnalysis
#' @references %% ~put references to the literature/web site here ~
#' @keywords misc
#' @examples
#'     findgRNAs(inputFilePath = system.file("extdata",
#'         "inputseq.fa", package = "CRISPRseek"),
#'         pairOutputFile = "testpairedgRNAs.xls",
#'         findPairedgRNAOnly = TRUE)               
#' 
#' 
#'     ##### predict gRNA efficacy using CRISPRscan
#'     featureWeightMatrixFile <- system.file("extdata", "Morenos-Mateo.csv",
#'             package = "CRISPRseek")
#' 
#'     findgRNAs(inputFilePath = system.file("extdata",
#'         "inputseq.fa", package = "CRISPRseek"),
#'         pairOutputFile = "testpairedgRNAs.xls",
#'         findPairedgRNAOnly = FALSE,
#'         calculategRNAEfficacy= TRUE,
#'         rule.set = "CRISPRscan", 
#'         baseBeforegRNA = 6, baseAfterPAM = 6,
#'         featureWeightMatrixFile = featureWeightMatrixFile,
#'         efficacyFile = "testCRISPRscanEfficacy.xls"
#'      )
#' 
#'      findgRNAs(inputFilePath = system.file("extdata",
#'         "testCRISPRscan.fa", package = "CRISPRseek"),
#'         pairOutputFile = "testpairedgRNAs.xls",
#'         findPairedgRNAOnly = FALSE,
#'         calculategRNAEfficacy= TRUE,
#'         rule.set = "CRISPRscan",
#'         baseBeforegRNA = 6, baseAfterPAM = 6,
#'         featureWeightMatrixFile = featureWeightMatrixFile,
#'         efficacyFile = "testCRISPRscanEfficacy.xls"
#'      )
#' 
#'     findgRNAs(inputFilePath = system.file("extdata", 
#'         "cpf1.fa", package = "CRISPRseek"), 
#'         findPairedgRNAOnly=FALSE, 
#'         pairOutputFile = "testpairedgRNAs-cpf1.xls", 
#'         PAM="TTTN", PAM.location = "5prime", PAM.size = 4, 
#'         overlap.gRNA.positions = c(19, 23),
#'         baseBeforegRNA = 8, baseAfterPAM = 26,
#'         calculategRNAEfficacy= TRUE, 
#'         rule.set = "DeepCpf1",
#'        efficacyFile = "testcpf1Efficacy.xls")
#' 
#'     findgRNAs(inputFilePath = system.file("extdata", 
#'              "cpf1.fa", package = "CRISPRseek"), 
#'              findPairedgRNAOnly=FALSE, 
#'              pairOutputFile = "testpairedgRNAs-cpf1.xls", 
#'              PAM="TTTN", PAM.location = "5prime", PAM.size = 4, 
#'              overlap.gRNA.positions = c(19,23),
#'              baseBeforegRNA = 8, baseAfterPAM = 26,
#'              calculategRNAEfficacy= TRUE, 
#'              rule.set = "DeepCpf1",
#'              efficacyFile = "testcpf1Efficacy.xls", baseEditing =  TRUE, 
#'              editingWindow=20, targetBase = "X")
#' 
#'     findgRNAs(inputFilePath = system.file("extdata", 
#'              "cpf1.fa", package = "CRISPRseek"), 
#'              findPairedgRNAOnly=FALSE, 
#'              pairOutputFile = "testpairedgRNAs-cpf1.xls", 
#'              PAM="TTTN", PAM.location = "5prime", PAM.size = 4, 
#'              overlap.gRNA.positions = c(19, 23),
#'              baseBeforegRNA = 8, baseAfterPAM = 26,
#'              calculategRNAEfficacy= TRUE, 
#'              rule.set = "DeepCpf1",
#'              efficacyFile = "testcpf1Efficacy.xls", baseEditing =  TRUE,
#'              editingWindow=20, targetBase = "C")
#' 
#'      inputSeq <-  DNAStringSet(paste(
#' "CCAGTTTGTGGATCCTGCTCTGGTGTCCTCCACACCAGAATCAGGGATCGAAAACTCA",
#' "TCAGTCGATGCGAGTCATCTAAATTCCGATCAATTTCACACTTTAAACG", sep =""))
#'      gRNAs <-  findgRNAs(inputFilePath =  inputSeq,
#'          pairOutputFile = "testpairedgRNAs1.xls",
#'          PAM.size = 3L,
#'          gRNA.size = 20L,
#'          overlap.gRNA.positions = c(17L,18L),
#'          PBS.length = 15,
#'          corrected.seq = "T",
#'          RT.template.pattern = "D$",
#'          RT.template.length = 8:30, 
#'          targeted.seq.length.change = 0,
#'          bp.after.target.end = 15,
#'          target.start = 46,
#'          target.end = 46,
#'          paired.orientation = "PAMin", min.gap = 20, max.gap = 90, 
#'          primeEditing = TRUE, findPairedgRNAOnly = TRUE)
#'
#' @importFrom Biostrings readDNAStringSet reverseComplement DNAStringSet
#' @importFrom BiocGenerics do.call rbind lapply unlist table cbind 
#' @importFrom IRanges Views IRanges findOverlaps width
#' @importFrom seqinr s2c
#' @importFrom S4Vectors orderIntegerPairs queryHits subjectHits 
#' @importFrom utils write.table read.csv
#' @export 

findgRNAs <-
    function (inputFilePath, 
        baseEditing = FALSE, targetBase = "C", editingWindow = 4:8,
        format = "fasta", PAM = "NGG", PAM.size = 3,
        findPairedgRNAOnly = FALSE, annotatePaired = TRUE,
        paired.orientation = c("PAMout", "PAMin"),
        enable.multicore = FALSE, n.cores.max = 6,
        gRNA.pattern = "", gRNA.size = 20, 
	overlap.gRNA.positions = c(17,18),
primeEditing = FALSE, 
        PBS.length = 13L,
        RT.template.length = 8:28,
        RT.template.pattern = "D$",
        corrected.seq,
        targeted.seq.length.change,
        bp.after.target.end = 15L,
        target.start,
        target.end,
primeEditingPaired.output = "pairedgRNAsForPE.xls",
        min.gap = 0, max.gap = 20, pairOutputFile, name.prefix = "",
	featureWeightMatrixFile = system.file("extdata", 
           "DoenchNBT2014.csv", package = "CRISPRseek"), 
        baseBeforegRNA = 4, 
        baseAfterPAM = 3,
	calculategRNAEfficacy = FALSE, efficacyFile,
        PAM.location = "3prime",
        rule.set = c("Root_RuleSet1_2014", 
            "Root_RuleSet2_2016", "CRISPRscan", "DeepCpf1"),
        chrom_acc)
{
    rule.set <- match.arg(rule.set)
    paired.orientation <- match.arg(paired.orientation)
    if (missing(inputFilePath)) {
        stop("inputFilePath containing the searching sequence 
	   or a DNAStringSet object is required!")
    }
    cut.site <- min(overlap.gRNA.positions)
    if (class(inputFilePath) != "DNAStringSet")
    {
    	if (! file.exists(inputFilePath)) {
           stop("inputfile specified as ", inputFilePath, " does not exists!")
        }
        if (format == "fasta" || format == "fastq")
        {
            subjects <- readDNAStringSet(inputFilePath, format, use.names = TRUE)
        }
        else
        {
            stop("format needs to be either fasta or fastq !")
        }
    }
    else
    {
	    subjects <- inputFilePath
    }
    PAM <- translatePattern(PAM)
    PAM <- paste("(?=", PAM, ")", sep="")

    gRNA.pattern <- translatePattern(gRNA.pattern)
    min.subject <- gRNA.size + PAM.size 
    subjects <-  subjects[width(subjects) >= min.subject,]
    if (length(subjects) == 0)
    {
        stop("The input file contains no sequence! This could be caused by 
            wrong format of the file. If file is created in mac, you could 
            reformat to text by typing tr \"\\r\" \"\\n\" >newfile in the 
            command line")
    }
    toAppend = FALSE
    colNames = TRUE
    names(subjects) <- gsub("'", "", names(subjects))
    names(subjects) <- gsub(" ", "", names(subjects))
    names(subjects) <- gsub("\t", ":", names(subjects))
     
    all.gRNAs.df <- do.call(rbind, lapply(1:length(subjects), function(p)
   {
        subject <- subjects[[p]]
        subjectname <- names(subjects)[p]
        revsubject <- reverseComplement(subject)
        plus.gRNAs <- 
           .getgRNA.cut.sites(subject, subjectname, PAM = PAM, 
               gRNA.pattern = gRNA.pattern, 
               gRNA.size = gRNA.size,
               cut.site = cut.site,
               PAM.size = PAM.size, 
               calculategRNAEfficacy = calculategRNAEfficacy,
               baseBeforegRNA = baseBeforegRNA,
               baseAfterPAM = baseAfterPAM,
               PAM.location = PAM.location,
               baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow)
       
       minus.gRNAs <-
           .getgRNA.cut.sites(revsubject, subjectname, PAM = PAM,
               gRNA.pattern = gRNA.pattern,
               gRNA.size = gRNA.size,
               cut.site = cut.site,
               PAM.size = PAM.size,
               calculategRNAEfficacy = calculategRNAEfficacy,
               baseBeforegRNA = baseBeforegRNA,
               baseAfterPAM = baseAfterPAM,
               reverse.subject = TRUE,
               PAM.location = PAM.location,
               baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow)

        if (length(plus.gRNAs) > 1)
            n.plus.gRNAs <- dim(plus.gRNAs)[1]
        else
            n.plus.gRNAs <- 0
        if (length(minus.gRNAs) > 1)
            n.minus.gRNAs <- dim(minus.gRNAs)[1]
        else
            n.minus.gRNAs <- 0
        if (n.minus.gRNAs == 0 && n.plus.gRNAs == 0)
        {
            cat(paste("No gRNAs found in the input sequence", subjectname))
            all.gRNAs <- DNAStringSet()
        }
        else if (findPairedgRNAOnly && (n.minus.gRNAs  * n.plus.gRNAs) == 0)
        {
            cat(paste("No paired gRNAs found in the input sequence", subjectname))
            all.gRNAs <- DNAStringSet()
        }
        else if (annotatePaired || findPairedgRNAOnly)
	    {
            temp = matrix(nrow =0, ncol=5)
            colnames(temp)[1:5] <- c( "ReversegRNAPlusPAM",
                "ReversegRNAName", "ForwardgRNAPlusPAM",
                "ForwardgRNAName", "gap")
	        if (p >1)
	        {
	            toAppend = TRUE
	            colNames = FALSE
	        }
            write.table(temp, file = pairOutputFile, sep = "\t",
                row.names = FALSE, quote = FALSE, append = toAppend, 
	            col.names = colNames)
            if (n.minus.gRNAs > 0 && n.plus.gRNAs > 0)
            {
                plus_start <- as.numeric(as.character(plus.gRNAs[ , 3L]))
                minus_end <- as.numeric(as.character(minus.gRNAs[ , 3L]))
                plus_end <- plus_start + gRNA.size + PAM.size - 1L
                minus_start <- minus_end - gRNA.size - PAM.size + 1L
		if (paired.orientation == "PAMout")
		{
                      pair_index <- .compute_pair_index3(plus_start,
                                                   minus_end,
                                                   min.gap, max.gap)
                      plus.index <- pair_index[[1L]]
                      minus.index <- pair_index[[2L]]
            	      gap <- pair_index[[3]]
                }
		else
		{
		     pair_index <- .compute_pair_index3(minus_start,
                                                   plus_end,
                                                   min.gap, max.gap)
                     plus.index <- pair_index[[2L]]
                     minus.index <- pair_index[[1L]]
                     gap <- pair_index[[3]]
		}
#cat("plus.index:", plus.index, "minus.index:", minus.index, "gap:", gap);
#cat("minus_start:", minus_start, "plus_start:", plus_start, "minus_end:", minus_end, "plus_end:", plus_end);
                if (!findPairedgRNAOnly)
                {
                    all.gRNAs <- DNAStringSet(c(plus.gRNAs[,1], minus.gRNAs[,1]))
                    names(all.gRNAs) <- c(plus.gRNAs[,2], minus.gRNAs[,2])
                    forEffi <- rbind(plus.gRNAs, minus.gRNAs)
                }
                if (length(minus.index) > 0 && length(plus.index) > 0)
                {
                    paired <- cbind(minus.gRNAs[minus.index,1], 
                        minus.gRNAs[minus.index, 2], plus.gRNAs[plus.index,1], 
                        plus.gRNAs[plus.index, 2],
                        gap = gap)
                    colnames(paired)[1:5] <- c( "ReversegRNAPlusPAM", 
                        "ReversegRNAName", "ForwardgRNAPlusPAM", 
                        "ForwardgRNAName", "gap")
                    if(findPairedgRNAOnly)
                    {
                        plus.index <- unique(plus.index)
                        minus.index <- unique(minus.index)
                        all.gRNAs <- DNAStringSet(c(plus.gRNAs[plus.index, 1], 
                            minus.gRNAs[minus.index,1]))
                        names(all.gRNAs) <- c(plus.gRNAs[plus.index,2], 
                            minus.gRNAs[minus.index,2])
                        forEffi <- rbind(plus.gRNAs[plus.index,], minus.gRNAs[minus.index,])
                    }
                    if (primeEditing)
                    {
                        paired <- designPEs(subject,
                           PAM.size = PAM.size,
                           gRNA.size = gRNA.size,
                           overlap.gRNA.positions = overlap.gRNA.positions,
                           PBS.length = PBS.length,
                           paired.gRNAs = data.frame(paired),
                           RT.template.length = RT.template.length,
                           RT.template.pattern = RT.template.pattern,
                           corrected.seq = corrected.seq,
                           targeted.seq.length.change = targeted.seq.length.change,
                           bp.after.target.end = bp.after.target.end,
                           target.start = target.start,
                           target.end = target.end,
                           primeEditingPaired.output =  primeEditingPaired.output, 
                           col.names = colNames, append = toAppend)
                         paired <- paired[,1:5]
                         all.gRNAs <- DNAStringSet(c(as.character(paired[,3]),
                            as.character(paired[,1])))
                         names(all.gRNAs) <- c(as.character(paired[,4]),
                            as.character(paired[,2])) 
                         all.gRNAs <- unique(all.gRNAs)
                         forEffi <- forEffi[forEffi[,2] %in% names(all.gRNAs),] 
                    }
                    if (dim(paired)[1] == 1)
                    {
                        write.table(paired, file = pairOutputFile, sep = "\t", 
                            row.names = FALSE, quote = FALSE, append = toAppend,
			                col.names = colNames)
                     }
                    else
                    {
                        write.table(paired[order(as.character(paired[,4])), ], 
                            file = pairOutputFile, sep = "\t", row.names = FALSE,
                            quote = FALSE, append = toAppend, col.names = colNames)
                    }
                } ## if paired found
                else if (findPairedgRNAOnly)
                {
                    cat(paste("No paired gRNAs found for sequence",
                        subjectname))
                   all.gRNAs <- DNAStringSet()
                   forEffi <- ""
                }
            }### if plus.gRNAs and minus.gRNAs not empty
            else if (!findPairedgRNAOnly && n.plus.gRNAs > 0)
            {
                all.gRNAs <- DNAStringSet(plus.gRNAs[,1])
                names(all.gRNAs) <- plus.gRNAs[,2]
                forEffi <- plus.gRNAs
            }
            else if (!findPairedgRNAOnly && n.minus.gRNAs > 0)
            {
                all.gRNAs <- DNAStringSet(minus.gRNAs[,1])
                names(all.gRNAs) <- minus.gRNAs[,2]
                forEffi <- minus.gRNAs
            }
        } ### annotatePaired and (no paired found or findPairedOnly)
        else if (n.minus.gRNAs > 0 && n.plus.gRNAs > 0)
        {
             all.gRNAs <- DNAStringSet(c(plus.gRNAs[,1], minus.gRNAs[,1]))
             names(all.gRNAs) <- c(plus.gRNAs[,2], minus.gRNAs[,2])
             forEffi <- rbind(plus.gRNAs, minus.gRNAs)
        }
        else if (n.minus.gRNAs > 0 )
        {
             all.gRNAs <- DNAStringSet(minus.gRNAs[,1])
             names(all.gRNAs) <- minus.gRNAs[,2]
             forEffi <- minus.gRNAs
        }
        else if (n.plus.gRNAs > 0 )
        {
            all.gRNAs <- DNAStringSet(plus.gRNAs[,1])
            names(all.gRNAs) <- plus.gRNAs[,2]
            forEffi <- plus.gRNAs
        }
        else
        {
            all.gRNAs <- DNAStringSet()
        }
        #if (length(all.gRNAs) == 0 && !findPairedgRNAOnly)
        #  cat(paste("No gRNAs found in the input sequence",
        #       subjectname))
        if (length(all.gRNAs) >0)
	       forEffi
    })) ## do call subjects
    if (calculategRNAEfficacy && length(all.gRNAs.df) > 4)
    {
        featureWeightMatrix <- read.csv(featureWeightMatrixFile, header=TRUE)
        if(rule.set == "Root_RuleSet1_2014")
        {
          effi <- calculategRNAEfficiency(all.gRNAs.df[,5], 
	    baseBeforegRNA = baseBeforegRNA, 
	    featureWeightMatrix = featureWeightMatrix, 
            enable.multicore = enable.multicore,
            n.cores.max = n.cores.max,
            gRNA.size = gRNA.size) 
        }
        else if (rule.set == "Root_RuleSet2_2016")
        {
          effi <- calculategRNAEfficiency2(all.gRNAs.df[,5])
        }
        else if (rule.set == "CRISPRscan")
        {
          effi <- calculategRNAEfficiencyCRISPRscan(all.gRNAs.df[,5], 
              featureWeightMatrix = featureWeightMatrix)
        }
        else if (rule.set == "DeepCpf1")
        {
           if (missing(chrom_acc))
                effi <- round(deepCpf1(extendedSequence = all.gRNAs.df[,5]), 3)
           else
		effi <- round(deepCpf1(extendedSequence = all.gRNAs.df[,5], 
                 chrom_acc = chrom_acc), 3)
        }
        extendedSequences <- cbind(all.gRNAs.df, effi)
        colnames(extendedSequences)  <- c("gRNAplusPAM", "name", "start", "strand", 
	    "extendedSequence", "gRNAefficacy")
        if (PAM.location == "3prime")
            extendedSequences[nchar(extendedSequences[,5]) 
	        < baseBeforegRNA + gRNA.size + PAM.size + baseAfterPAM, 6] <- 
                "extended sequence too short"
        else
            extendedSequences[nchar(extendedSequences[,5])
                < baseBeforegRNA + baseAfterPAM , 6] <-
                "extended sequence too short"
        write.table(extendedSequences,
             file = efficacyFile, sep="\t", row.names = FALSE)
    }
    #else
    #{
#	extendedSequences <- all.gRNAs.df[,1:4]
#	colnames(extendedSequences)  = c("gRNAplusPAM", "name", "start", "strand")
#    }
    all.gRNAs <- DNAStringSet(all.gRNAs.df[,1])
    names(all.gRNAs) = all.gRNAs.df[,2]
    #list(all.gRNAs=all.gRNAs, extendedSequences = extendedSequences)
    all.gRNAs
}
