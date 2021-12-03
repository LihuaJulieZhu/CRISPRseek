#' Compare 2 input sequences/sequence sets for possible guide RNAs (gRNAs)
#' 
#' Generate all possible guide RNAs (gRNAs) for two input sequences, or two
#' sets of sequences, and generate scores for potential off-targets in the other
#' sequence.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param inputFile1Path Sequence input file 1 path that contains one of the
#' two sequences to be searched for potential gRNAs. It can also be a
#' DNAStringSet object with names field set. Please see examples below.
#' @param inputFile2Path Sequence input file 2 path that contains one of the
#' two sequences to be searched for potential gRNAs. It can also be a
#' DNAStringSet object with names field set. Please see examples below.
#' @param inputNames Name of the input sequences when inputFile1Path and
#' inputFile2Path are DNAStringSet instead of file path
#' @param format Format of the input files, fasta, fastq and bed format are
#' supported, default fasta
#' @param header Indicate whether the input file contains header, default
#' FALSE, only applies to bed format
#' @param findgRNAsWithREcutOnly Indicate whether to find gRNAs overlap with
#' restriction enzyme recognition pattern
#' @param searchDirection Indicate whether perfrom gRNA in both sequences and
#' off-target search against each other (both) or search gRNA in input1 and
#' off-target analysis in input2 (1to2), or vice versa (2to1)
#' @param BSgenomeName BSgenome object. Please refer to available.genomes in
#' BSgenome package. For example, BSgenome.Hsapiens.UCSC.hg19 for hg19,
#' BSgenome.Mmusculus.UCSC.mm10 for mm10, BSgenome.Celegans.UCSC.ce6 for ce6,
#' BSgenome.Rnorvegicus.UCSC.rn5 for rn5, BSgenome.Drerio.UCSC.danRer7 for Zv9,
#' and BSgenome.Dmelanogaster.UCSC.dm3 for dm3
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
#' @param editingWindow.offtargets Applicable only when baseEditing is set to
#' TRUE. It is used to indicate the effective editing window to consider for
#' the offtargets search only, default to 4 to 8 which is for the original CBE
#' system. Please change it accordingly if the system you use have a different
#' editing window, or you would like to include offtargets with the target base
#' in a larger editing window.
#' @param REpatternFile File path containing restriction enzyme cut patters
#' @param minREpatternSize Minimum restriction enzyme recognition pattern
#' length required for the enzyme pattern to be searched for, default 6
#' @param findgRNAs Indicate whether to find gRNAs from the sequences in the
#' input file or skip the step of finding gRNAs, default TRUE for both input
#' sequences.  Set it to FALSE if the input file contains user selected gRNAs
#' plus PAM already.
#' @param removegRNADetails Indicate whether to remove the detailed gRNA
#' information such as efficacy file and restriction enzyme cut sites, default
#' false for both input sequences. Set it to TRUE if the input file contains
#' the user selected gRNAs plus PAM already.
#' @param exportAllgRNAs Indicate whether to output all potential gRNAs to a
#' file in fasta format, genbank format or both. Default to no.
#' @param annotatePaired Indicate whether to output paired information, default
#' to FALSE
#' @param overlap.gRNA.positions The required overlap positions of gRNA and
#' restriction enzyme cut site, default 17 and 18
#' @param findPairedgRNAOnly Choose whether to only search for paired gRNAs in
#' such an orientation that the first one is on minus strand called reverse
#' gRNA and the second one is on plus strand called forward gRNA. TRUE or
#' FALSE, default FALSE
#' @param min.gap Minimum distance between two oppositely oriented gRNAs to be
#' valid paired gRNAs. Default 0
#' @param max.gap Maximum distance between two oppositely oriented gRNAs to be
#' valid paired gRNAs. Default 20
#' @param gRNA.name.prefix The prefix used when assign name to found gRNAs,
#' default _gR, short for guided RNA.
#' @param PAM.size PAM length, default 3
#' @param gRNA.size The size of the gRNA, default 20
#' @param PAM PAM sequence after the gRNA, default NGG
#' @param PAM.pattern Regular expression of PAM, default NNG or NGN for spCas9.
#' For cpf1, ^TTTN since it is a 5 prime PAM sequence
#' @param allowed.mismatch.PAM Maximum number of mismatches allowed to the PAM
#' sequence, default to 1 for PAM.pattern NNG or NGN PAM
#' @param max.mismatch Maximum mismatch allowed to search the off targets in
#' the other sequence, default 3
#' @param outputDir the directory where the sequence comparison results will be
#' written to
#' @param upstream upstream offset from the bed input starts to search for gRNA
#' and/or offtargets, default 0
#' @param downstream downstream offset from the bed input ends to search for
#' gRNA and/or offtargets, default 0
#' @param weights numeric vector size of gRNA length, default c(0, 0, 0.014, 0,
#' 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828,
#' 0.615, 0.804, 0.685, 0.583) which is used in Hsu et al., 2013 cited in the
#' reference section
#' @param overwrite overwrite the existing files in the output directory or
#' not, default TRUE
#' @param baseBeforegRNA Number of bases before gRNA used for calculating gRNA
#' efficiency, default 4 Please note, for PAM located on the 5 prime, need to
#' specify the number of bases before the PAM sequence plus PAM size.
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
#' @param foldgRNAs Default FALSE. If set to TRUE, summary file will contain
#' minimum free energy of the secondary structure of gRNA with gRNA backbone
#' from GeneRfold package provided that GeneRfold package has been installed.
#' @param gRNA.backbone gRNA backbone constant region sequence. Default to the
#' sequence in Sp gRNA backbone.
#' @param temperature temperature in celsius. Default to 37 celsius.
#' @param scoring.method Indicates which method to use for offtarget cleavage
#' rate estimation, currently two methods are supported, Hsu-Zhang and CFDscore
#' @param subPAM.activity Applicable only when scoring.method is set to
#' CFDscore A hash to represent the cleavage rate for each alternative sub PAM
#' sequence relative to preferred PAM sequence
#' @param subPAM.position Applicable only when scoring.method is set to
#' CFDscore The start and end positions of the sub PAM. Default to 22 and 23
#' for SP with 20bp gRNA and NGG as preferred PAM
#' @param PAM.location PAM location relative to gRNA. For example, spCas9 PAM
#' is located on the 3 prime (3prime) while cpf1 PAM is located on the 5 prime
#' (5prime)
#' @param rule.set Specify a rule set scoring system for calculating gRNA
#' efficacy. Please note that Root_RuleSet2_2016 requires the following python
#' packages with specified verion and python 2.7.  1. scikit-learn 0.16.1 2.
#' pickle 3. pandas 4. numpy 5. scipy
#' @param mismatch.activity.file Applicable only when scoring.method is set to
#' CFDscore A comma separated (csv) file containing the cleavage rates for all
#' possible types of single nucleotide mismatche at each position of the gRNA.
#' By default, using the supplemental Table 19 from Doench et al., Nature
#' Biotechnology 2016
#' @return Return a data frame with all potential gRNAs from both sequences. In
#' addition, a tab delimited file scoresFor2InputSequences.xls is also saved in
#' the outputDir, sorted by scoreDiff descending. 
#' \item{name}{name of the gRNA}
#' \item{gRNAPlusPAM}{gRNA plus PAM sequence}
#' \item{targetInSeq1}{target/off-target sequence including PAM in the 1st
#' input sequence file} \item{targetInSeq2}{target/off-target sequence incuding
#' PAM in the 2nd input sequence file}
#' \item{guideAlignment2Offtarget}{alignment of gRNA to the other input
#' sequence (off-target sequence)} \item{offTargetStrand}{strand of the other
#' sequence (off-target sequence) the gRNA align to} \item{scoreForSeq1}{score
#' for the target sequence in the 1st input sequence file}
#' \item{scoreForSeq2}{score for the target sequence in the 1st input sequence
#' file} \item{mismatch.distance2PAM}{distances of mismatch to PAM, e.g., 14
#' means the mismatch is 14 bp away from PAM} \item{n.mismatch}{number of
#' mismatches between the off-target and the gRNA} \item{targetSeqName}{the
#' name of the input sequence where the target sequence is located}
#' \item{scoreDiff}{scoreForSeq1 - scoreForSeq2} \item{bracket.notation}{folded
#' gRNA in bracket notation} \item{mfe.sgRNA}{minimum free energy of sgRNA}
#' \item{mfe.diff}{mfe.sgRNA-mfe.backbone} \item{mfe.backbone}{minimum free
#' energy of the gRNA backbone by itself}
#' @note %% ~~further notes~~
#' @author Lihua Julie Zhu
#' @seealso CRISPRseek
#' @references Patrick D Hsu, David A Scott, Joshua A Weinstein, F Ann Ran,
#' Silvana Konermann, Vineeta Agarwala, Yinqing Li, Eli J Fine, Xuebing Wu,
#' Ophir Shalem, Thomas J Cradick, Luciano A Marraffini, Gang Bao & Feng Zhang
#' (2013) DNA targeting specificity of rNA-guided Cas9 nucleases. Nature
#' Biotechnology 31:827-834
#' @keywords misc
#' @examples
#' 
#'     library(CRISPRseek)
#'     inputFile1Path <- system.file("extdata", "rs362331T.fa",
#'             package = "CRISPRseek")
#'     inputFile2Path <- system.file("extdata", "rs362331C.fa",
#'             package = "CRISPRseek")
#'     REpatternFile <- system.file("extdata", "NEBenzymes.fa", 
#'             package = "CRISPRseek")
#'     seqs <- compare2Sequences(inputFile1Path, inputFile2Path,
#'         outputDir = getwd(), 
#'         REpatternFile = REpatternFile, overwrite = TRUE)
#' 
#'     seqs2 <- compare2Sequences(inputFile1Path, inputFile2Path,
#'                inputNames=c("Seq1", "Seq2"),
#'                scoring.method = "CFDscore",
#'                outputDir = getwd(), 
#'                overwrite = TRUE, baseEditing = TRUE)
#' 
#'     inputFile1Path <- 
#' DNAStringSet(
#' "TAATATTTTAAAATCGGTGACGTGGGCCCAAAACGAGTGCAGTTCCAAAGGCACCCACCTGTGGCAG"
#' )
#'     ## when set inputFile1Path to a DNAStringSet object, it is important
#'     ## to call names
#'     names(inputFile1Path) <- "seq1"
#'     
#'     inputFile2Path <- 
#' DNAStringSet(
#' "TAATATTTTAAAATCGGTGACGTGGGCCCAAAACGAGTGCAGTTCCAAAGGCACCCACCTGTGGCAG"
#' )
#'      ## when set inputFile2Path to a DNAStringSet object, it is important 
#'     ## to call names
#' 
#'     names(inputFile2Path) <- "seq2"
#' 
#'     seqs <- compare2Sequences(inputFile1Path, inputFile2Path,
#'           inputNames=c("Seq1", "Seq2"),
#'           scoring.method = "CFDscore",
#'           outputDir = getwd(), 
#'           overwrite = TRUE)
#' 
#'     seqs2 <- compare2Sequences(inputFile1Path, inputFile2Path,
#'                inputNames=c("Seq1", "Seq2"),
#'                scoring.method = "CFDscore",
#'                outputDir = getwd(), 
#'                overwrite = TRUE, baseEditing = TRUE)
#' @importFrom GenomicRanges intersect
#' @importFrom BiocGenerics rbind
#' @importFrom Biostrings readDNAStringSet DNAStringSet
#' @importFrom utils read.csv read.table write.table 
#' @importFrom hash hash
#' @export
compare2Sequences <- function(inputFile1Path, inputFile2Path, inputNames=c("Seq1", "Seq2"), 
	format = c("fasta", "fasta"), header = FALSE, findgRNAsWithREcutOnly = FALSE,
    searchDirection = c("both","1to2", "2to1"), BSgenomeName,
    baseEditing = FALSE, targetBase = "C", editingWindow = 4:8,
    editingWindow.offtargets = 4:8,
    REpatternFile=system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek"),
    minREpatternSize = 6, findgRNAs = c(TRUE, TRUE), removegRNADetails = c(FALSE, FALSE), 
    exportAllgRNAs = c("no", "all", "fasta", "genbank"), annotatePaired =  FALSE,
    overlap.gRNA.positions = c(17, 18), findPairedgRNAOnly = FALSE, 
    min.gap = 0, max.gap = 20, gRNA.name.prefix = "_gR", PAM.size = 3, 
    gRNA.size = 20, PAM = "NGG", PAM.pattern = "NNG$|NGN$",
    allowed.mismatch.PAM = 1, max.mismatch = 3, 
    outputDir, upstream =0, downstream = 0,
    weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 
    0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 
    0.685, 0.583), overwrite = FALSE, baseBeforegRNA = 4, 
    baseAfterPAM = 3, featureWeightMatrixFile = system.file("extdata", 
       "DoenchNBT2014.csv", package = "CRISPRseek"), foldgRNAs = FALSE, 
        gRNA.backbone="GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU",
        temperature = 37,
        scoring.method = c("Hsu-Zhang", "CFDscore"),
        subPAM.activity = hash( AA =0,
          AC =   0,
          AG = 0.259259259,
          AT = 0,
          CA = 0,
          CC = 0,
          CG = 0.107142857,
          CT = 0,
          GA = 0.069444444,
          GC = 0.022222222,
          GG = 1,
          GT = 0.016129032,
          TA = 0,
          TC = 0,
          TG = 0.038961039,
          TT = 0),
     subPAM.position = c(22, 23),
     PAM.location = "3prime",
     rule.set = c("Root_RuleSet1_2014", "Root_RuleSet2_2016", "CRISPRscan", "DeepCpf1"),
     mismatch.activity.file = system.file("extdata", 
         "NatureBiot2016SuppTable19DoenchRoot.csv", 
         package = "CRISPRseek")
    )
{
        scoring.method <- match.arg(scoring.method)
        exportAllgRNAs <- match.arg(exportAllgRNAs)
        searchDirection <- match.arg(searchDirection)
        rule.set <- match.arg(rule.set)
        PAM.p.letters <- strsplit(PAM.pattern, split="")[[1]]
        if (PAM.location == "3prime" && PAM.p.letters[length(PAM.p.letters)] != "$")
           PAM.pattern <- paste0(PAM.pattern, "$")
        if (PAM.location == "5prime" && PAM.p.letters[1] != "^")
           PAM.pattern <- paste0("^", PAM.pattern)
        if (scoring.method == "Hsu-Zhang")
        {
             if (length(weights) !=  gRNA.size)
                 stop("Please make sure the size of weights vector 
                     equals to the gRNA.size!\n")
        }
        else if (scoring.method ==  "CFDscore")
        {
            mismatch.activity <- read.csv(mismatch.activity.file)
            required.col <- c("Mismatch.Type", "Position", "Percent.Active")
            if (length(intersect(colnames(mismatch.activity), required.col)) !=
                length(required.col))
                stop("Please rename the mismatch activity file column to contain at least
                   these 3 column names: Mismatch.Type, Position, Percent.Active\n")
        }
        if (class(inputFile1Path) != "DNAStringSet" || class(inputFile2Path) != "DNAStringSet")
        {
	    if ((format[1] == "bed" || format[2] == "bed") && 
                (missing(BSgenomeName) || class(BSgenomeName) != "BSgenome"))
                   stop("BSgenomeName is required as BSgenome object when input file is in bed format!")
        }
	append = ifelse(overwrite, FALSE, TRUE)
	if (class(inputFile1Path) != "DNAStringSet")
	{
	    outputDir1 <- file.path(outputDir, paste(basename(inputFile1Path),
               format(Sys.time(), "%b-%d-%Y"), sep="-"))
	}
	else
	{
           outputDir1 <- file.path(outputDir, "File1")
	}
	if (class(inputFile2Path) != "DNAStringSet")
	{
	   outputDir2 <- file.path(outputDir, paste(basename(inputFile2Path),
               format(Sys.time(), "%b-%d-%Y"), sep="-"))
	}
	else
	{
	   outputDir2 <- file.path(outputDir, "File2")
	}
	if(searchDirection == "both" || searchDirection == "1to2")
	{
		cat("search for gRNAs for input file1...\n")
                if (findgRNAs[1])
                {
		    tryCatch(
			 (gRNAs1 = offTargetAnalysis(inputFile1Path, format = format[1], 
			     findgRNAs = findgRNAs[1], annotatePaired =  annotatePaired,
                             exportAllgRNAs = exportAllgRNAs, gRNAoutputName= inputNames[1],
			     findPairedgRNAOnly = findPairedgRNAOnly, chromToSearch = "", 
			     findgRNAsWithREcutOnly = findgRNAsWithREcutOnly, 
			     REpatternFile = REpatternFile, minREpatternSize = minREpatternSize, 
			     overlap.gRNA.positions =  overlap.gRNA.positions,
			     min.gap = min.gap, max.gap = max.gap, 
			     gRNA.name.prefix = gRNA.name.prefix, PAM.size = PAM.size,
			     gRNA.size = gRNA.size, PAM = PAM, PAM.pattern = PAM.pattern,
                             outputUniqueREs = FALSE,
                             baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow,
                             editingWindow.offtargets = editingWindow.offtargets,
			     outputDir = outputDir1, upstream.search = upstream,
			     downstream.search = downstream,
			     weights = weights, foldgRNAs = FALSE, overwrite = overwrite,
			     featureWeightMatrixFile = featureWeightMatrixFile, 
            		     baseBeforegRNA = baseBeforegRNA, BSgenomeName = BSgenomeName,
            		     baseAfterPAM = baseAfterPAM, header = header,
                             subPAM.position = subPAM.position,
                             subPAM.activity = subPAM.activity,
                             PAM.location = PAM.location,
                             rule.set = rule.set,
                             mismatch.activity.file = mismatch.activity.file)), 
			 error = function(e) {print(e); gRNAs1 = DNAStringSet()})
             }
             else if (class(inputFile1Path) == "DNAStringSet") 
             {
                gRNAs1 <- inputFile1Path
             }
            else  if (format[1] == "fasta" || format[1] == "fastq")
            {
                gRNAs1 <- readDNAStringSet(inputFile1Path, format[1],
                      use.names = TRUE)
            }
            else
            {
                stop("format needs to be either fasta or fastq for gRNA file!")
            }
	}
	if(searchDirection == "both" || searchDirection == "2to1")
	{
             if(findgRNAs[2])
             {
		cat("search for gRNAs for input file2...\n")
		tryCatch((gRNAs2 = offTargetAnalysis(inputFile2Path, format = format[2],          
			findgRNAs = findgRNAs[2], annotatePaired =  annotatePaired,
                        exportAllgRNAs = exportAllgRNAs, gRNAoutputName = inputNames[2],
			findPairedgRNAOnly = findPairedgRNAOnly, chromToSearch = "",
            		findgRNAsWithREcutOnly = findgRNAsWithREcutOnly, 
            		REpatternFile = REpatternFile, minREpatternSize = minREpatternSize,
                        outputUniqueREs = FALSE,
            		overlap.gRNA.positions =  overlap.gRNA.positions, 
            		min.gap = min.gap, max.gap = max.gap, 
            		gRNA.name.prefix = gRNA.name.prefix, PAM.size = PAM.size,
            		gRNA.size = gRNA.size, PAM = PAM, PAM.pattern = PAM.pattern, 
                        baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow,
                        editingWindow.offtargets = editingWindow.offtargets,
            		outputDir = outputDir2, upstream.search = upstream,
			downstream.search = downstream, 
            		weights = weights, foldgRNAs = FALSE, overwrite = overwrite,
                        featureWeightMatrixFile = featureWeightMatrixFile,
                        baseBeforegRNA = baseBeforegRNA, BSgenomeName = BSgenomeName,
                        baseAfterPAM = baseAfterPAM, header = header, 
                        subPAM.position = subPAM.position,
                        subPAM.activity = subPAM.activity,
                        PAM.location = PAM.location,
                        rule.set = rule.set,
                        mismatch.activity.file = mismatch.activity.file)), 
			error=function(e) {print(e); gRNAs2 = DNAStringSet()})
            }
            else if (class(inputFile2Path) == "DNAStringSet")
            {
                gRNAs2 <- inputFile2Path
            }
            else  if (format[2] == "fasta" || format[2] == "fastq")
            {
                gRNAs2 <- readDNAStringSet(inputFile2Path, format[2],
                      use.names = TRUE)
            }
            else
            {
                stop("format needs to be either fasta or fastq for gRNA file!")
            }
	}
    print("Scoring ...")
    if (class(inputFile1Path) != "DNAStringSet")
    {
	if (format[1] == "bed")
	    subjects1 <- getSeqFromBed(inputFile1Path, header = header,
                BSgenomeName = BSgenomeName, upstream = upstream, downstream = downstream)
	else
    	    subjects1 <- readDNAStringSet(inputFile1Path, format = format[1],
    	        use.names = TRUE)
    }
    else
    {
	   subjects1 <- inputFile1Path
    }
   if (class(inputFile2Path) != "DNAStringSet")
   {
	if (format[2] == "bed")
            subjects2 <- getSeqFromBed(inputFile2Path, header = header, 
	        BSgenomeName = BSgenomeName, upstream = upstream, downstream = downstream)
	else
    	    subjects2 <- readDNAStringSet(inputFile2Path, format = format[2],
       	        use.names = TRUE)
   }
   else
   {
	   subjects2 <- inputFile2Path
   }
    outfile <- tempfile(tmpdir = getwd())
    seqname <- names(subjects2)
    seqname <- gsub("'", "", seqname)
    seqname <- gsub(" ", "", seqname)
    seqname <- gsub("\t", ":", seqname)
    names(subjects2) <- seqname
    #write.table(subjects2, file="testSeqs")
    #revsubject <- reverseComplement(subjects2)
    #chrom.len <- nchar(as.character(subjects2))
    if(searchDirection == "both" || searchDirection == "1to2")
	{
	   for (j in 1:length(subjects2))
	   {
                if (j == 1)
                    hits <- searchHits(gRNAs = gRNAs1, PAM = PAM, PAM.pattern = PAM.pattern, 
                        seqs = subjects2[[j]], seqname = names(subjects2)[j],
                        max.mismatch = max.mismatch, PAM.size = PAM.size, 
                        gRNA.size = gRNA.size, allowed.mismatch.PAM = allowed.mismatch.PAM,
                        PAM.location = PAM.location, outfile = outfile, 
                        baseEditing = baseEditing, targetBase = targetBase,
                        editingWindow = editingWindow.offtargets) 
                else
                    hits <- rbind(hits, searchHits(gRNAs = gRNAs1, PAM = PAM, 
                        PAM.pattern = PAM.pattern,
                        seqs = subjects2[[j]], seqname = names(subjects2)[j],
                        max.mismatch = max.mismatch, PAM.size = PAM.size,
                        gRNA.size = gRNA.size, allowed.mismatch.PAM = allowed.mismatch.PAM,
                        PAM.location = PAM.location, outfile = outfile,
                        baseEditing = baseEditing, targetBase = targetBase,
                        editingWindow = editingWindow.offtargets))
           }
	} # end of if searchDirection == "both" or searchDirection == "1to2"
	cat("finish off-target search in sequence 2\n") 
   	seqname <- names(subjects1)
	seqname <- gsub("'", "", seqname)
	seqname <- gsub(" ", "", seqname)
	seqname <- gsub("\t", ":", seqname)
	#revsubject <- reverseComplement(subjects1)
	#revsubject <- reverseComplement(subjects1[[1]])
	names(subjects1) <- seqname
	#chrom.len = nchar(as.character(subjects1))
	if(searchDirection == "both" || searchDirection == "2to1")
	{
	    for (j in 1:length(subjects1))
	    {
                if (j == 1  && !exists("hits"))
                    hits <- searchHits(gRNAs = gRNAs2, PAM = PAM, PAM.pattern = PAM.pattern, 
                       seqs = subjects1[[j]], seqname = names(subjects1)[j],
                       max.mismatch = max.mismatch, PAM.size = PAM.size, 
                       gRNA.size = gRNA.size, allowed.mismatch.PAM = allowed.mismatch.PAM,
                       PAM.location = PAM.location, outfile = outfile,
                       baseEditing = baseEditing, targetBase = targetBase,
                       editingWindow = editingWindow.offtargets) 
               else
                   hits <- rbind(hits, searchHits(gRNAs = gRNAs2, PAM = PAM, 
                       PAM.pattern = PAM.pattern,
                       seqs = subjects1[[j]], seqname = names(subjects1)[j],
                       max.mismatch = max.mismatch, PAM.size = PAM.size, 
                       gRNA.size = gRNA.size, allowed.mismatch.PAM = allowed.mismatch.PAM,
                       PAM.location = PAM.location, outfile = outfile,
                       baseEditing = baseEditing, targetBase = targetBase,
                       editingWindow = editingWindow.offtargets))
	    }
	} # if searchDirection == "both" or searchDirection == "2to1"
	cat("finish off-target search in sequence 1\n")
        if (removegRNADetails[1])
            unlink(outputDir1, recursive = TRUE)
        if (removegRNADetails[2])
            unlink(outputDir2, recursive = TRUE)
	featureVectors <- buildFeatureVectorForScoring(hits = hits, 
	    canonical.PAM = PAM, gRNA.size = gRNA.size,
            subPAM.position = subPAM.position, PAM.location = PAM.location)
	cat("finish feature vector building\n")
        if ( scoring.method ==  "CFDscore")
            scores <- getOfftargetScore2(featureVectors,
                subPAM.activity = subPAM.activity,
                mismatch.activity.file = mismatch.activity.file)
        else
	    scores <- getOfftargetScore(featureVectors, weights = weights)
	cat("finish score calculation\n")
	targetInSeq1 <- scores$gRNAPlusPAM
	targetInSeq2 <- scores$gRNAPlusPAM
        if ( scoring.method ==  "CFDscore")
        {	
            scoreForSeq1 <- rep(1, dim(scores)[1])
	    scoreForSeq2 <- rep(1, dim(scores)[1])
            max.score <- 1
        }
        else
        {
            scoreForSeq1 <- rep(100, dim(scores)[1])
            scoreForSeq2 <- rep(100, dim(scores)[1])
            max.score <- 100
        }
	if(length(subjects1) == 1 && length(subjects2) == 1)
	{
		scoreForSeq1[scores$chrom == names(subjects1)] <- 
			scores$score[scores$chrom == names(subjects1)]
		scoreForSeq2[scores$chrom == names(subjects2)] <- 
			scores$score[scores$chrom == names(subjects2)]
		targetInSeq1[scores$chrom == names(subjects1)] <- 
			scores$OffTargetSequence[scores$chrom == names(subjects1)]
		targetInSeq2[scores$chrom == names(subjects2)] <- 
			scores$OffTargetSequence[scores$chrom == names(subjects2)]
		seqNames <- c(names(subjects1), names(subjects2))
		targetSeqName = unlist(lapply(
			scores$chrom, function(i) {seqNames[seqNames !=i]}))
	}
	else if (searchDirection == "1to2")
	{
		scoreForSeq2 = scores$score
		targetInSeq2 = scores$OffTargetSequence
		targetSeqName = scores$name
	}
	else if (searchDirection == "2to1")
	{
		scoreForSeq1 = scores$score
		targetInSeq1 = scores$OffTargetSequence
		targetSeqName = scores$name
	}
	else
	{	
		seqNames <- c(names(subjects1), names(subjects2))
		fileIndex <- c(rep(1, length(subjects1)), rep(2, length(subjects2)))
		offTargetFiles <- fileIndex[match(scores$chrom,seqNames)]
		targetSeqName <- scores$name
		scoreForSeq1[offTargetFiles == 1] <- 
			scores[offTargetFiles == 1, ]$score
		scoreForSeq2[offTargetFiles == 2] <- 
			scores[offTargetFiles == 2, ]$score
		targetInSeq1[offTargetFiles ==1] <- 
			scores[offTargetFiles == 1,]$OffTargetSequence
		targetInSeq2[offTargetFiles == 2] <- 
			scores[offTargetFiles == 2,]$OffTargetSequence
	}
	targetSeqName <- gsub(paste(gRNA.name.prefix, "[0-9]+[f|r]", sep=""), "", targetSeqName)
	targetSeqName <- gsub(paste("_gR", "[0-9]+[f|r]", sep=""), "", targetSeqName)
	targetSeqName <- gsub( "Start[0-9]+End[0-9]+", "", targetSeqName)		
	seqs.new <- as.data.frame(cbind(name = scores$name,
		gRNAPlusPAM = scores$gRNAPlusPAM,
		targetInSeq1 = targetInSeq1,
		targetInSeq2 = targetInSeq2,
	    guideAlignment2OffTarget = scores$alignment,
	    offTargetStrand = scores$strand,
		scoreForSeq1 = scoreForSeq1,
		scoreForSeq2 = scoreForSeq2,
		mismatch.distance2PAM = as.character(scores$mismatch.distance2PAM),
		n.mismatch = scores$n.mismatch,
		offTarget = scores$forViewInUCSC,
		targetSeqName = targetSeqName)
	)
	 
	if(searchDirection == "both" || searchDirection == "1to2")
	{
	if (length(setdiff(names(gRNAs1),seqs.new[,1])) >0)
	{
		gRNAnames <- names(gRNAs1)[!names(gRNAs1) %in% seqs.new[,1]]
		targetSeqName <- gsub(paste(gRNA.name.prefix, "[0-9]+[f|r]", sep=""), "", gRNAnames)
		targetSeqName <- gsub(paste("_gR", "[0-9]+[f|r]", sep=""), "", targetSeqName)
		targetSeqName <- gsub( "Start[0-9]+End[0-9]+", "", targetSeqName)		
		seqs1.only <- as.data.frame(cbind(name = gRNAnames,
			gRNAPlusPAM = as.character(gRNAs1)[names(gRNAs1) %in% gRNAnames],
			targetInSeq1 = as.character(gRNAs1)[names(gRNAs1) %in% gRNAnames],
			targetInSeq2 = rep("NA", length(targetSeqName)),
		  guideAlignment2OffTarget = rep("NA", length(targetSeqName)),
		  offTargetStrand = rep("NA", length(targetSeqName)),
			scoreForSeq1 = rep(max.score, length(targetSeqName)),
			scoreForSeq2 = rep(0, length(targetSeqName)),
			mismatch.distance2PAM = rep("NA", length(targetSeqName)),
			n.mismatch = rep("NA", length(targetSeqName)),
			offTarget = rep("NA", length(targetSeqName)),
			targetSeqName = targetSeqName)
			)
		seqs.new <- rbind(seqs.new, seqs1.only)
	}
	}
	if(searchDirection == "both" || searchDirection == "2to1")
	{
	if (length(setdiff(names(gRNAs2),seqs.new[,1])) >0)
	{
		gRNAnames <- names(gRNAs2)[!names(gRNAs2) %in% seqs.new[,1]]
		targetSeqName <- gsub(paste(gRNA.name.prefix, "[0-9]+[f|r]", sep=""), "",gRNAnames)
		targetSeqName <- gsub(paste("_gR", "[0-9]+[f|r]", sep=""), "", targetSeqName)
		targetSeqName <- gsub( "_Start[0-9]+End[0-9]+", "", targetSeqName)
		seqs2.only <- as.data.frame(cbind(name = gRNAnames,
    		gRNAPlusPAM = as.character(gRNAs2)[names(gRNAs2) %in% gRNAnames],
		    targetInSeq1 = rep("NA", length(targetSeqName)),
			  targetInSeq2 = as.character(gRNAs2)[names(gRNAs2) %in% gRNAnames],
		    guideAlignment2OffTarget = rep("NA", length(targetSeqName)),
        offTargetStrand = rep("NA", length(targetSeqName)),
        scoreForSeq1 = rep(0, length(targetSeqName)),
			  scoreForSeq2 = rep(max.score, length(targetSeqName)),
			  mismatch.distance2PAM = rep("NA", length(targetSeqName)),
			  n.mismatch = rep("NA", length(targetSeqName)),
			  offTarget = rep("NA", length(targetSeqName)),
			  targetSeqName = targetSeqName)
		)
		seqs.new <- rbind(seqs.new, seqs2.only)
	}
	}
	seqs <- unique(cbind(seqs.new, gRNAefficacy = 0, 
        scoreDiff = round(as.numeric(as.character(seqs.new[,7])) - 
                            as.numeric(as.character(seqs.new[,8]),4), 3)))
#	rownames(seqs) <- seqs[,1]
#	seqs <- as.data.frame(seqs)
#save(seqs, file="seqs")
        if (substr(outputDir1, nchar(outputDir1), nchar(outputDir1)) != .Platform$file.sep)
    	{
       	    outputDir1 <- paste(outputDir1, "", sep = .Platform$file.sep)
    	}
	if (substr(outputDir2, nchar(outputDir2), nchar(outputDir2)) != .Platform$file.sep)
        {
            outputDir2 <- paste(outputDir2, "", sep = .Platform$file.sep)
        }
	if (searchDirection == "both")
	{
            if (findgRNAs[1])
            {
		eff1File <- paste(outputDir1, "gRNAefficacy.xls", sep = "")
                gRNAeff1 <- read.table(eff1File,sep="\t", header=TRUE, stringsAsFactors=FALSE)
            }  
            if (findgRNAs[2])
            {
		eff2File <- paste(outputDir2, "gRNAefficacy.xls", sep = "")
		gRNAeff2 <- read.table(eff2File,sep="\t", header=TRUE, stringsAsFactors=FALSE)
            }
            if (findgRNAs[1])
            {
                if (findgRNAs[2])
		    gRNAeff <- rbind(gRNAeff1, gRNAeff2) 
                else
                    gRNAeff <- gRNAeff1
            }
            else if (findgRNAs[2])
                gRNAeff <- gRNAeff2
        }
	if(searchDirection == "1to2" && findgRNAs[1])
	{
		eff1File <- paste(outputDir1, "gRNAefficacy.xls", sep = "")
		gRNAeff <- read.table(eff1File,sep="\t", header=TRUE, stringsAsFactors=FALSE)
        }
	if(searchDirection == "2to1" && findgRNAs[2])
	{
		eff2File <- paste(outputDir2, "gRNAefficacy.xls", sep = "")
		gRNAeff <- read.table(eff2File,sep="\t", header=TRUE, stringsAsFactors=FALSE)
	}
        if(exists("gRNAeff"))
        {
	    m <- match(seqs$name, gRNAeff$name)
	    seqs$gRNAefficacy <- gRNAeff$gRNAefficacy[m]
        }
	originalDir <- getwd()
	setwd(outputDir)
        if (foldgRNAs)
        {
           gRNAs.withoutPAM <- substr(as.character(seqs$gRNAPlusPAM), 1, gRNA.size)
           folded.gRNAs <- foldgRNAs(gRNAs.withoutPAM, gRNA.backbone = gRNA.backbone, 
           temperature = temperature)
           if (length(dim(folded.gRNAs)) > 0)
               seqs <- cbind(seqs, folded.gRNAs[,-1])
        }
	if (dim(seqs)[1] ==1)
	{
		write.table(seqs, file = "scoresFor2InputSequences.xls",
			sep = "\t", row.names = FALSE, col.names=TRUE)
	}
	else
	{
             write.table(seqs[order(as.numeric(seqs[,dim(seqs)[2]]), decreasing = TRUE), ], 
	        file = "scoresFor2InputSequences.xls",
		sep = "\t", row.names = FALSE, col.names=TRUE)
	}
	print("Done!")
        if (removegRNADetails[1])
            unlink(outputDir1, recursive = TRUE)
        if (removegRNADetails[2])
            unlink(outputDir2, recursive = TRUE)
	setwd(originalDir)
#scores
	seqs
}
