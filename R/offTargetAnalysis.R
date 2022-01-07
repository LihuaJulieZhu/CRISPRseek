#' Design target-specific guide RNAs for CRISPR-Cas9 system in one function
#'
#' Design target-specific guide RNAs (gRNAs) and predict relative indel
#' fequencies for CRISPR-Cas9 system by automatically calling findgRNAs,
#' filtergRNAs, searchHits, buildFeatureVectorForScoring, getOfftargetScore,
#' filterOfftarget, calculating gRNA cleavage efficiency, and predict gRNA
#' efficacy, indels and their frequencies.
#'
#' %% ~~ If necessary, more details than the description above ~~
#'
#' @param inputFilePath Sequence input file path or a DNAStringSet object that
#' contains sequences to be searched for potential gRNAs
#' @param format Format of the input file, fasta, fastq and bed are supported,
#' default fasta
#' @param header Indicate whether the input file contains header, default
#' FALSE, only applies to bed format
#' @param gRNAoutputName Specify the name of the gRNA outupt file when
#' inputFilePath is DNAStringSet object instead of file path
#' @param findgRNAs Indicate whether to find gRNAs from the sequences in the
#' input file or skip the step of finding gRNAs, default TRUE. Set it to FALSE
#' if the input file contains user selected gRNAs plus PAM already.
#' @param exportAllgRNAs Indicate whether to output all potential gRNAs to a
#' file in fasta format, genbank format or both. Default to both.
#' @param findgRNAsWithREcutOnly Indicate whether to find gRNAs overlap with
#' restriction enzyme recognition pattern
#' @param REpatternFile File path containing restriction enzyme cut patterns
#' @param minREpatternSize Minimum restriction enzyme recognition pattern
#' length required for the enzyme pattern to be searched for, default 4
#' @param overlap.gRNA.positions The required overlap positions of gRNA and
#' restriction enzyme cut site, default 17 and 18. For Cpf1, you can set it to 19 and 23.
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
#' @param min.gap Minimum distance between two oppositely oriented gRNAs to be
#' valid paired gRNAs. Default 0
#' @param enable.multicore Indicate whether enable parallel processing, default
#' FALSE. For super long sequences with lots of gRNAs, suggest set it to TRUE
#' @param n.cores.max Indicating maximum number of cores to use in multi core
#' mode, i.e., parallel processing, default 6. Please set it to 1 to disable
#' multicore processing for small dataset.
#' @param max.gap Maximum distance between two oppositely oriented gRNAs to be
#' valid paired gRNAs. Default 20
#' @param gRNA.name.prefix The prefix used when assign name to found gRNAs,
#' default gRNA, short for guided RNA.
#' @param PAM.size PAM length, default 3
#' @param gRNA.size The size of the gRNA, default 20
#' @param PAM PAM sequence after the gRNA, default NGG
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
#' @param chromToExclude Specify the chromosome not to search. If specified as
#' "", meaning to search chromosomes specified by chromToSearch. By default, to
#' exclude haplotype blocks from offtarget search in hg19, i.e., chromToExclude
#' = c("chr17_ctg5_hap1","chr4_ctg9_hap1", "chr6_apd_hap1", "chr6_cox_hap2",
#' "chr6_dbb_hap3", "chr6_mann_hap4", "chr6_mcf_hap5","chr6_qbl_hap6",
#' "chr6_ssto_hap7")
#' @param max.mismatch Maximum mismatch allowed in off target search, default
#' 3. Warning: will be considerably slower if set > 3
#' @param PAM.pattern Regular expression of protospacer-adjacent motif (PAM),
#' default NNG$|NGN$ for spCas9. For cpf1, ^TTTN since it is a 5 prime PAM
#' sequence
#' @param allowed.mismatch.PAM Maximum number of mismatches allowed in the PAM
#' sequence for offtarget search, default to 1 to allow NGN and NNG PAM pattern
#' for offtarget identification.
#' @param gRNA.pattern Regular expression or IUPAC Extended Genetic Alphabet to
#' represent gRNA pattern, default is no restriction. To specify that the gRNA
#' must start with GG for example, then set it to ^GG. Please see
#' help(translatePattern) for a list of IUPAC Extended Genetic Alphabet.
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
#' the offtargets search only, default to 4 to 8 (1 means the most distal site
#' from the 3' PAM, the most proximla site from the 5' PAM), which is for the
#' original CBE system.  Please change it accordingly if the system you use
#' have a different editing window, or you would like to include offtargets
#' with the target base in a larger editing window.
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
#' @param min.score minimum score of an off target to included in the final
#' output, default 0
#' @param topN top N off targets to be included in the final output, default
#' 1000
#' @param topN.OfftargetTotalScore top N off target used to calculate the total
#' off target score, default 10
#' @param annotateExon Choose whether or not to indicate whether the off target
#' is inside an exon or not, default TRUE
#' @param txdb TxDb object, for creating and using TxDb object, please refer to
#' GenomicFeatures package. For a list of existing TxDb object, please search
#' for annotation package starting with Txdb at
#' http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData,
#' such as
#' \itemize{
#' \item{TxDb.Rnorvegicus.UCSC.rn5.refGene} - for rat
#' \item{TxDb.Mmusculus.UCSC.mm10.knownGene} - for mouse
#' \item{TxDb.Hsapiens.UCSC.hg19.knownGene} - for human
#' \item{TxDb.Dmelanogaster.UCSC.dm3.ensGene} - for Drosophila
#' \item{TxDb.Celegans.UCSC.ce6.ensGene} - for C.elegans
#' }
#' @param orgAnn organism annotation mapping such as org.Hs.egSYMBOL in
#' org.Hs.eg.db package for human
#' @param ignore.strand default to TRUE when annotating to gene
#' @param outputDir the directory where the off target analysis and reports
#' will be written to
#' @param fetchSequence Fetch flank sequence of off target or not, default TRUE
#' @param upstream upstream offset from the off target start, default 200
#' @param downstream downstream offset from the off target end, default 200
#' @param upstream.search upstream offset from the bed input starts to search
#' for gRNAs, default 0
#' @param downstream.search downstream offset from the bed input ends to search
#' for gRNAs, default 0
#' @param weights Applicable only when scoring.method is set to Hsu-Zhang a
#' numeric vector size of gRNA length, default c(0, 0, 0.014, 0, 0, 0.395,
#' 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615,
#' 0.804, 0.685, 0.583) which is used in Hsu et al., 2013 cited in the
#' reference section
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
#' @param useScore Default TRUE, display in gray scale with the darkness
#' indicating the gRNA efficacy.  The taller bar shows the Cas9 cutting site.
#' If set to False, efficacy will not show.  Instead, gRNAs in plus strand will
#' be colored red and gRNAs in negative strand will be colored green.
#' @param useEfficacyFromInputSeq Default FALSE. If set to TRUE, summary file
#' will contain gRNA efficacy calculated from input sequences instead of from
#' off-target analysis. Set it to TRUE if the input sequence is from a
#' different species than the one used for off-target analysis.
#' @param outputUniqueREs Default TRUE. If set to TRUE, summary file will
#' contain REs unique to the cleavage site within 100 or 200 bases surrounding
#' the gRNA sequence.
#' @param foldgRNAs Default FALSE. If set to TRUE, summary file will contain
#' minimum free energy of the secondary structure of gRNA with gRNA backbone
#' from GeneRfold package provided that GeneRfold package has been installed.
#' @param gRNA.backbone gRNA backbone constant region sequence. Default to the
#' sequence in Sp gRNA backbone.
#' @param temperature temperature in celsius. Default to 37 celsius.
#' @param overwrite overwrite the existing files in the output directory or
#' not, default FALSE
#' @param scoring.method Indicates which method to use for offtarget cleavage
#' rate estimation, currently two methods are supported, Hsu-Zhang and CFDscore
#' @param subPAM.activity Applicable only when scoring.method is set to
#' CFDscore A hash to represent the cleavage rate for each alternative sub PAM
#' sequence relative to preferred PAM sequence
#' @param subPAM.position Applicable only when scoring.method is set to
#' CFDscore The start and end positions of the sub PAM. Default to 22 and 23
#' for spCas9 with 20bp gRNA and NGG as preferred PAM. For Cpf1, it could be
#' c(1,2).
#' @param PAM.location PAM location relative to gRNA. For example, default to
#' 3prime for spCas9 PAM.  Please set to 5prime for cpf1 PAM since it's PAM is
#' located on the 5 prime end
#' @param rule.set Specify a rule set scoring system for calculating gRNA
#' efficacy. Please note that Root_RuleSet2_2016 requires the following python
#' packages with specified verion and python 2.7.  1. scikit-learn 0.16.1 2.
#' pickle 3. pandas 4. numpy 5. scipy
#' @param chrom_acc Optional binary variable indicating chromatin accessibility
#' information with 1 indicating accessible and 0 not accessible.
#' @param calculategRNAefficacyForOfftargets Default to TRUE to output gRNA
#' efficacy for offtargets as well as ontargets. Set it to FALSE if only need
#' gRNA efficacy calculated for ontargets only to speed up the analysis. Please
#' refer to https://support.bioconductor.org/p/133538/#133661 for potential use
#' cases of offtarget efficacies.
#' @param mismatch.activity.file Applicable only when scoring.method is set to
#' CFDscore A comma separated (csv) file containing the cleavage rates for all
#' possible types of single nucleotide mismatche at each position of the gRNA.
#' By default, using the supplemental Table 19 from Doench et al., Nature
#' Biotechnology 2016
#' @param predIndelFreq Default to FALSE. Set it to TRUE to output the
#' predicted indels and their frequencies.
#' @param predictIndelFreq.onTargetOnly Default to TRUE, indicating that indels
#' and their frequencies will be predicted for ontargets only. Usually,
#' researchers are only interested in predicting the editing outcome for the
#' ontargets since any editing in the offtargets are unwanted. Set it to FALSE
#' if you are interested in predicting indels and their frequencies for
#' offtargets. It will take longer time to run if you set it to FALSE.
#' @param method.indelFreq Currently only Lindel method has been implemented.
#' Please let us know if you think additional methods should be made available.
#' Lindel is compatible with both Python2.7 and Python3.5 or higher. Please
#' type help(predictRelativeFreqIndels) to get more details.
#' @param baseBeforegRNA.indelFreq Default to 13 for Lindel method.
#' @param baseAfterPAM.indelFreq Default to 24 for Lindel method.
#' @return Four tab delimited files are generated in the output directory:
#' \item{OfftargetAnalysis.xls}{ - detailed information of off targets}
#' \item{Summary.xls}{ - summary of the gRNAs}
#' \item{REcutDetails.xls}{ - restriction enzyme cut sites of each gRNA}
#' \item{pairedgRNAs.xls}{ - potential paired gRNAs}
#' @note %% ~~further notes~~
#' @author Lihua Julie Zhu
#' @seealso CRISPRseek
#' @references Patrick D Hsu, David A Scott, Joshua A Weinstein, F Ann Ran,
#' Silvana Konermann, Vineeta Agarwala, Yinqing Li, Eli J Fine, Xuebing Wu,
#' Ophir Shalem, Thomas J Cradick, Luciano A Marraffini, Gang Bao & Feng Zhang
#' (2013) DNA targeting specificity of rNA-guided Cas9 nucleases. Nature
#' Biotechnology 31:827-834
#'
#' Doench JG, Hartenian E, Graham DB, Tothova Z, Hegde M, Smith I, Sullender M,
#' Ebert BL, Xavier RJ, Root DE. Rational design of highly active sgRNAs for
#' CRISPR-Cas9-mediated gene inactivation. Nat Biotechnol. 2014 Sep 3. doi:
#' 10.1038 nbt.3026
#'
#' Lihua Julie Zhu, Benjamin R. Holmes, Neil Aronin and Michael Brodsky.
#' CRISPRseek: a Bioconductor package to identify target-specific guide RNAs
#' for CRISPR-Cas9 genome-editing systems. Plos One Sept 23rd 2014
#'
#' Moreno-Mateos, M., Vejnar, C., Beaudoin, J. et al. CRISPRscan: designing
#' highly efficient sgRNAs for CRISPR-Cas9 targeting in vivo.  Nat Methods 12,
#' 982–988 (2015) doi:10.1038/nmeth.3543
#'
#' Doench JG et al., Optimized sgRNA design to maximize activity and minimize
#' off-target effects of CRISPR-Cas9. Nature Biotechnology Jan 18th 2016
#'
#' Anzalone et al., Search-and-replace genome editing without double-strand
#' breaks or donor DNA. Nature October 2019
#' https://www.nature.com/articles/s41586-019-1711-4
#'
#' Wei Chen, Aaron McKenna, Jacob Schreiber et al., Massively parallel
#' profiling and predictive modeling of the outcomes of CRISPR/Cas9-mediated
#' double-strand break repair, Nucleic Acids Research, Volume 47, Issue 15, 05
#' September 2019, Pages 7989–8003, https://doi.org/10.1093/nar/gkz487
#'
#' Kim et al., Deep learning improves prediction of CRISPR–Cpf1
#' guide RNA activityNat Biotechnol 36, 239–241 (2018).
#' https://doi.org/10.1038/nbt.4061
#'
#' @keywords misc
#' @examples
#'
#' 	library(CRISPRseek)
#' 	library("BSgenome.Hsapiens.UCSC.hg19")
#' 	library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' 	library(org.Hs.eg.db)
#' 	outputDir <- getwd()
#' 	inputFilePath <- system.file("extdata", "inputseq.fa",
#'             package = "CRISPRseek")
#' 	REpatternFile <- system.file("extdata", "NEBenzymes.fa",
#'             package = "CRISPRseek")
#' 	results <- offTargetAnalysis(inputFilePath, findgRNAsWithREcutOnly = TRUE,
#'             REpatternFile = REpatternFile, findPairedgRNAOnly = FALSE,
#'             annotatePaired = FALSE,
#'             BSgenomeName = Hsapiens, chromToSearch = "chrX",
#'             txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#' 	    orgAnn = org.Hs.egSYMBOL, max.mismatch = 1,
#'             outputDir = outputDir, overwrite = TRUE)
#'
#'        #### predict indels and their frequecies for target sites
#'
#'        if (interactive())
#'        {
#'           results <- offTargetAnalysis(inputFilePath,findgRNAsWithREcutOnly = TRUE,
#'             findPairedgRNAOnly = FALSE,
#'             annotatePaired = FALSE,
#'             BSgenomeName = Hsapiens, chromToSearch = "chrX",
#'             txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
#' 	    orgAnn = org.Hs.egSYMBOL, max.mismatch = 1,
#'             outputDir = outputDir, overwrite = TRUE,
#'             predIndelFreq=TRUE, predictIndelFreq.onTargetOnly= TRUE)
#'
#'           names(results$indelFreq)
#'           head(results$indelFreq[[1]])
#'           ### save the indel frequences to tab delimited files, one file for each target/offtarget site.
#'           mapply(write.table, results$indelFreq, file=paste0(names(results$indelFreq), '.xls'), sep = "\t", row.names = FALSE)
#'
#'        #### predict gRNA efficacy using CRISPRscan
#'        featureWeightMatrixFile <- system.file("extdata", "Morenos-Mateo.csv",
#'             package = "CRISPRseek")
#'
#'        results <- offTargetAnalysis(inputFilePath, findgRNAsWithREcutOnly = TRUE,
#'             REpatternFile = REpatternFile, findPairedgRNAOnly = FALSE,
#'             annotatePaired = FALSE,
#'             BSgenomeName = Hsapiens, chromToSearch = "chrX",
#'             txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'             orgAnn = org.Hs.egSYMBOL, max.mismatch = 1,
#'             rule.set = "CRISPRscan",
#'             baseBeforegRNA = 6, baseAfterPAM = 6,
#'             featureWeightMatrixFile = featureWeightMatrixFile,
#'             outputDir = outputDir, overwrite = TRUE)
#'
#'        ######## PAM is on the 5 prime side, e.g., Cpf1
#'        results <- offTargetAnalysis(inputFilePath = system.file("extdata",
#'               "cpf1-2.fa", package = "CRISPRseek"), findgRNAsWithREcutOnly =  FALSE,
#'           findPairedgRNAOnly = FALSE,
#'           annotatePaired = FALSE,
#'           BSgenomeName = Hsapiens,
#'           chromToSearch = "chr8",
#'           txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'           orgAnn = org.Hs.egSYMBOL, max.mismatch = 4,
#'           baseBeforegRNA = 8, baseAfterPAM = 26,
#'           rule.set = "DeepCpf1",
#'           overlap.gRNA.positions = c(19, 23),
#'           useEfficacyFromInputSeq = FALSE,
#'           outputDir = getwd(),
#'           overwrite = TRUE, PAM.location = "5prime",PAM.size = 4,
#'           PAM = "TTTN", PAM.pattern = "^TNNN", allowed.mismatch.PAM = 2,
#'           subPAM.position = c(1,2))
#'
#'         results1 <- offTargetAnalysis(inputFilePath, findgRNAsWithREcutOnly =  FALSE,
#'                  REpatternFile = REpatternFile, findPairedgRNAOnly = FALSE,
#'                  annotatePaired = FALSE,
#'                  BSgenomeName = Hsapiens, chromToSearch = "chrX",
#'                  txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                  orgAnn = org.Hs.egSYMBOL, max.mismatch = 4,
#'                  outputDir = outputDir, overwrite = TRUE, PAM.location = "5prime",
#'                  PAM = "TGT", PAM.pattern = "^T[A|G]N", allowed.mismatch.PAM = 2,
#'                  subPAM.position = c(1,2), baseEditing = TRUE, editingWindow =20, targetBase = "G")
#'
#'         results.testBE <- offTargetAnalysis(inputFilePath, findgRNAsWithREcutOnly =  FALSE,
#'                  REpatternFile = REpatternFile, findPairedgRNAOnly = FALSE,
#'                  annotatePaired = FALSE,
#'                  BSgenomeName = Hsapiens, chromToSearch = "chrX",
#'                  txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                  orgAnn = org.Hs.egSYMBOL, max.mismatch = 4,
#'                  outputDir = outputDir, overwrite = TRUE, PAM.location = "5prime",
#'                  PAM = "TGT", PAM.pattern = "^T[A|G]N", allowed.mismatch.PAM = 2,
#'                  subPAM.position = c(1,2), baseEditing = TRUE,
#'                  editingWindow = 10:20, targetBase = "A")
#'
#'         inputFilePath <- DNAStringSet(paste(
#' "CCAGTTTGTGGATCCTGCTCTGGTGTCCTCCACACCAGAATCAGGGATCGAAAA",
#' "CTCATCAGTCGATGCGAGTCATCTAAATTCCGATCAATTTCACACTTTAAACG", sep =""))
#'         names(inputFilePath) <- "testPE"
#'         results3 <- offTargetAnalysis(inputFilePath,
#'             gRNAoutputName = "testPEgRNAs",
#'             BSgenomeName = Hsapiens, chromToSearch = "chrX",
#'             txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'             orgAnn = org.Hs.egSYMBOL, max.mismatch = 1,
#'             outputDir = outputDir, overwrite = TRUE,
#'             PAM.size = 3L,
#'             gRNA.size = 20L,
#'             overlap.gRNA.positions = c(17L,18L),
#'             PBS.length = 15,
#'             corrected.seq = "T",
#'             RT.template.pattern = "D$",
#'             RT.template.length = 8:30,
#'             targeted.seq.length.change = 0,
#'             bp.after.target.end = 15,
#'             target.start = 20,
#'             target.end = 20,
#'             paired.orientation = "PAMin", min.gap = 20, max.gap = 90,
#'             primeEditing = TRUE, findPairedgRNAOnly = TRUE)
#'        }
#' @importFrom hash hash
#' @importFrom utils read.csv write.table read.table
#' @importFrom GenomicRanges intersect setdiff
#' @importFrom Biostrings writeXStringSet readDNAStringSet
#' @importFrom BiocGenerics rbind as.data.frame cbind unlist lapply
#' @export
offTargetAnalysis <-
    function(inputFilePath, format = "fasta", header=FALSE,
        gRNAoutputName, findgRNAs = TRUE,
        exportAllgRNAs = c("all", "fasta", "genbank", "no"),
        findgRNAsWithREcutOnly = FALSE,
	REpatternFile = system.file("extdata", "NEBenzymes.fa",
            package = "CRISPRseek"),
	minREpatternSize = 4,
	overlap.gRNA.positions = c(17, 18), findPairedgRNAOnly = FALSE,
        annotatePaired = TRUE, paired.orientation = c("PAMout","PAMin"),
        enable.multicore = FALSE, n.cores.max = 6,
        min.gap = 0, max.gap = 20, gRNA.name.prefix = "",
	PAM.size = 3, gRNA.size = 20, PAM = "NGG", BSgenomeName,
        chromToSearch = "all",
        chromToExclude = c("chr17_ctg5_hap1","chr4_ctg9_hap1", "chr6_apd_hap1",
"chr6_cox_hap2", "chr6_dbb_hap3", "chr6_mann_hap4", "chr6_mcf_hap5","chr6_qbl_hap6",
"chr6_ssto_hap7"),
	max.mismatch = 3,
        PAM.pattern = "NNG$|NGN$", allowed.mismatch.PAM = 1,
        gRNA.pattern = "",
        baseEditing = FALSE, targetBase = "C", editingWindow = 4:8,
        editingWindow.offtargets = 4:8,
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
        min.score = 0, topN = 1000,
        topN.OfftargetTotalScore = 10,
        annotateExon = TRUE, txdb, orgAnn, ignore.strand = TRUE, outputDir,
        fetchSequence = TRUE, upstream = 200, downstream = 200,
	upstream.search = 0, downstream.search = 0,
        weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445,
        0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583),
	baseBeforegRNA = 4, baseAfterPAM = 3,
	featureWeightMatrixFile = system.file("extdata", "DoenchNBT2014.csv",
		package = "CRISPRseek"),
	useScore = TRUE, useEfficacyFromInputSeq = FALSE,
	outputUniqueREs = TRUE, foldgRNAs = FALSE,
        gRNA.backbone="GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU",
        temperature = 37,
        overwrite = FALSE,
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
     chrom_acc,
     calculategRNAefficacyForOfftargets = TRUE,
     mismatch.activity.file = system.file("extdata",
         "NatureBiot2016SuppTable19DoenchRoot.csv",
         package = "CRISPRseek"),
     predIndelFreq = FALSE,
     predictIndelFreq.onTargetOnly = TRUE,
     method.indelFreq = "Lindel",
     baseBeforegRNA.indelFreq = 13,
     baseAfterPAM.indelFreq = 24
)
{
    cat("Validating input ...\n")
    scoring.method <- match.arg(scoring.method)
    exportAllgRNAs <- match.arg(exportAllgRNAs)
    rule.set <- match.arg(rule.set)
    PAM.p.letters <- strsplit(PAM.pattern, split="")[[1]]
    if (PAM.location == "3prime" && PAM.p.letters[length(PAM.p.letters)] != "$")
      PAM.pattern <- paste0(PAM.pattern, "$")
    if (PAM.location == "5prime" && PAM.p.letters[1] != "^")
      PAM.pattern <- paste0("^", PAM.pattern)
    if (rule.set == "DeepCpf1")
    {
        baseBeforegRNA <- 8
        baseAfterPAM <- 26
        if (scoring.method == "CFDscore" && subPAM.activity$TT < 1)
          subPAM.activity = hash( AA =0,
              AC = 0,
              AG = 0,
          AT = 0.1,
          CA = 0,
          CC = 0,
          CG = 0,
          CT = 0.05,
          GA = 0,
          GC = 0,
          GG = 0,
          GT = 0.05,
          TA = 0.2,
          TC = 0.1,
          TG = 0.1,
          TT = 1)
    }
    else if (rule.set %in% c("Root_RuleSet1_2014",
        "Root_RuleSet2_2016", "CRISPRscan"))
    {
        if (PAM.location == "3prime")
        {
            baseBeforegRNA <- 4
            baseAfterPAM <- 3
        }
        else
        {
            baseBeforegRNA <- 4 + PAM.size
            baseAfterPAM <- 3 + gRNA.size
        }
    }
    if (scoring.method ==  "CFDscore")
    {
        mismatch.activity <- read.csv(mismatch.activity.file)
        required.col <- c("Mismatch.Type", "Position", "Percent.Active")
        if (length(intersect(colnames(mismatch.activity), required.col)) !=
            length(required.col))
           stop("Please rename the mismatch activity file column to contain at least
              these 3 column names: Mismatch.Type, Position, Percent.Active\n")
    }
    else if (scoring.method == "Hsu-Zhang")
    {
         if (length(weights) !=  gRNA.size)
             stop("Please make sure the size of weights vector
                 equals to the gRNA.size!\n")
    }
    if(findgRNAsWithREcutOnly && findgRNAs && !file.exists(REpatternFile))
    {
        stop("Please specify an REpattern file as fasta file with
            restriction enzyme recognition sequences!")
    }
    if (missing(inputFilePath)) {
        stop("inputFilePath containing the searching sequence, coordinate or a DNAStringSet
             object is required!")
    }
    if (substr(outputDir, nchar(outputDir), nchar(outputDir)) != .Platform$file.sep)
    {
        outputDir <- paste(outputDir, "", sep = .Platform$file.sep)
    }
    if (file.exists(outputDir) && ! overwrite)
    {
        cat(outputDir, "exists already. Please type 1 if you want to
            overwrite the outputDir and 2 if you want to exit.", fill = TRUE)
	input <- readline()
	if(input != 1) { stop("Please change the outputDir!") }
    }
    if (!file.exists(outputDir))
    {
        dir.create(outputDir)
    }
    if (annotatePaired || findPairedgRNAOnly)
       pairOutputFile <- paste(outputDir, "pairedgRNAs.xls", sep = "")
    REcutDetailFile <- paste(outputDir, "REcutDetails.xls", sep = "")
    bedFile<- paste(outputDir, "gRNAsCRISPRseek.bed", sep = "")
    if (missing(gRNAoutputName) && class(inputFilePath) == "DNAStringSet")
	    stop("Please enter a name for the gRNA ouput file name when DNAStringSet instead of file path provided!")
    if (class(inputFilePath) != "DNAStringSet" && missing(gRNAoutputName))
	    gRNAoutputName = strsplit(basename(inputFilePath), split=".",
		    fixed=TRUE)[[1]][1]
    if (format =="bed")
    {
        if (missing(BSgenomeName) || class(BSgenomeName) != "BSgenome") {
            stop("BSgenomeName is required as BSgenome object when input file is in bed format!")
        }
        inputFilePath <- getSeqFromBed(inputFilePath, header = header, BSgenomeName = BSgenomeName)
        #### format for filtergRNAs
        format <- "fasta"
    }
    if (findgRNAs)
    {
        cat("Searching for gRNAs ...\n")
	efficacyFile <- paste(outputDir, "gRNAefficacy.xls", sep = "")
	if (chromToSearch == "" || useEfficacyFromInputSeq)
            potential.gRNAs <- findgRNAs(inputFilePath,
               overlap.gRNA.positions = overlap.gRNA.positions,
               baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow,
               primeEditing = primeEditing,
               findPairedgRNAOnly = findPairedgRNAOnly,
               annotatePaired = annotatePaired,
               paired.orientation = paired.orientation,
               pairOutputFile = pairOutputFile, PAM = PAM,
               PAM.location = PAM.location,
               gRNA.pattern = gRNA.pattern, PAM.size = PAM.size,
               gRNA.size = gRNA.size, min.gap = min.gap,
               max.gap = max.gap, name.prefix = gRNA.name.prefix,
               format = format, featureWeightMatrixFile = featureWeightMatrixFile,
               baseBeforegRNA = baseBeforegRNA,
	       baseAfterPAM = baseAfterPAM ,
    	       calculategRNAEfficacy = TRUE, efficacyFile = efficacyFile,
               rule.set = rule.set, chrom_acc = chrom_acc)
         else
	    potential.gRNAs <- findgRNAs(inputFilePath,
               overlap.gRNA.positions = overlap.gRNA.positions,
              baseEditing = baseEditing, targetBase = targetBase, editingWindow = editingWindow,
              primeEditing = primeEditing,
              PBS.length = PBS.length,
               RT.template.length = RT.template.length,
               RT.template.pattern = RT.template.pattern,
               corrected.seq = corrected.seq,
               targeted.seq.length.change = targeted.seq.length.change,
               bp.after.target.end = bp.after.target.end,
               target.start = target.start,
               target.end = target.end,
               primeEditingPaired.output =  primeEditingPaired.output,
               findPairedgRNAOnly = findPairedgRNAOnly,
               annotatePaired = annotatePaired,
               paired.orientation = paired.orientation,
               enable.multicore = enable.multicore,
               n.cores.max = n.cores.max,
               pairOutputFile = pairOutputFile, PAM = PAM,
	       gRNA.pattern = gRNA.pattern, PAM.size = PAM.size,
               PAM.location = PAM.location,
               gRNA.size = gRNA.size, min.gap = min.gap,
               max.gap = max.gap, name.prefix = gRNA.name.prefix,
               format = format,  rule.set = rule.set, chrom_acc = chrom_acc)
	if (length(potential.gRNAs) == 0)
        {
		return(cat("no gRNAs found!"))
        }

	if (length(potential.gRNAs) > 0 && (exportAllgRNAs == "fasta" || exportAllgRNAs == "all"))
	{
		writeXStringSet(potential.gRNAs, filepath= file.path(outputDir,
                     paste(gRNAoutputName,"allgRNAs.fa", sep="")))
	}
	if (length(potential.gRNAs) > 0 && (exportAllgRNAs == "genbank" || exportAllgRNAs == "all"))
	{
		if (class(inputFilePath) == "DNAStringSet")
			subjects <- inputFilePath
		else
			subjects <- readDNAStringSet(inputFilePath, format=format,
				use.names = TRUE)
                names(subjects) <- gsub( "\t", "", names(subjects))
                names(subjects) <- gsub( "\n", "", names(subjects))
                names(subjects) <- gsub( " ", "", names(subjects))
   	        locuses <- names(subjects)

		names.gRNA <- names(potential.gRNAs)
		for (i in 1:length(locuses))
		{
			thisLocus <- gsub("'", "", locuses[i])
        		thisLocus <- gsub(" ", "", thisLocus)
			thisSeq <- tolower(as.character(subjects[[i]]))
			n.bp <- nchar(thisSeq)
			temp <- strsplit(names.gRNA, split=paste(
	         		thisLocus,"_gR",sep=""))
			locus <- paste("LOCUS       ", thisLocus,
                                        "                     ", n.bp,
					" bp    dna     linear   UNK", sep="")
			definition <- paste("DEFINITION  CRISPRseek output for ",
				    gRNAoutputName, " sequence", sep = "")
			accession <- "ACCESSION   unknown"
			features <- "FEATURES             Location/Qualifiers"
			header = rbind(locus, definition, accession, features)
			found.gRNA <- 0
			for (j in 1:length(temp))
			{
				if (length(temp[[j]]) >1){
					found.gRNA <- found.gRNA + 1
					if (found.gRNA == 1)
					{
					    thisFile <- file.path(outputDir,
					      paste(thisLocus, "gbk", sep="."))
                            		    write(header, thisFile)
					}
					if  (length(grep("f", temp[[j]])) >0)
					{
						temp1 <-strsplit(temp[[j]], "f")
						isForward <- TRUE
					}
					else
					{
						temp1 <-strsplit(temp[[j]], "r")
						isForward <- FALSE
					}
					feature <- temp1[[2]][2]
					feature[is.na(feature)] <- ""
					location <- temp1[[2]][1]
					if (isForward)
					{
				       	    Start <- location
					    End <- as.numeric(Start) + max(overlap.gRNA.positions) -
						min(overlap.gRNA.positions)
					    write(paste("     misc_bind       ", Start, "..",
                                                End, sep = ""), append = TRUE, sep="\n",
                                                file = thisFile)
					     write(paste("                     /note=\"gRNAf",
						as.character(feature),
                                                "\"", sep = ""), append = TRUE, sep="\n", file = thisFile)
					}
					else
					{
                                            End <- location
                                            Start <- as.numeric(End) - max(overlap.gRNA.positions) +
						min(overlap.gRNA.positions)
					    write(paste("     misc_bind       complement(",
						    Start, "..", End, ")", sep = ""),
						    append = TRUE, sep="\n", file = thisFile)
					    write(paste("                     /note=\"gRNAr",
						feature,
                            			"\"", sep = ""), append = TRUE, sep="\n", file = thisFile)
					}
				}
			}
			if (found.gRNA > 0){
			    write("ORIGIN", append = TRUE, sep="\n", file = thisFile)
                    	    seq.lines <- floor(nchar(thisSeq) / 60) + 1
                            for (k in 1:seq.lines) {
                                line.start <- (k - 1) * 60 + 1
                        	line.end <- min(line.start + 59, nchar(thisSeq))
                        	n.leading.spaces <- 9 - nchar(line.start)
                        	leading.spaces <- paste(rep(" ", n.leading.spaces),
                            		collapse = "")
                        	seq.thisLine <- substr(thisSeq, line.start, line.end)
                        	len.thisLine <- nchar(seq.thisLine)
                        	n.seg <- floor(len.thisLine /10) + 1
                        	for (l in 1:n.seg) {
                            		seg.start <- (l -1) * 10 + 1
                            		seg.end <- min(seg.start + 9, len.thisLine)
                            		if (l == 1)
                                		seq.thisLine.formatted <- substr(seq.thisLine,
                                    			seg.start, seg.end)
                            		else
                                	seq.thisLine.formatted <- paste(
                                    		seq.thisLine.formatted,
                                    		substr(seq.thisLine, seg.start, seg.end),
                                    		sep = " ")
                             	}
                        	write(paste(leading.spaces, line.start, " ",
                            		seq.thisLine.formatted, sep = ""),
                            		append = TRUE, sep="\n", file = thisFile)
                    	}
				write("//", append = TRUE, sep="\n", file = thisFile)
			}
		    }
		}
	if (findPairedgRNAOnly && length(potential.gRNAs) >0)
	{
	    gRNAs.RE <- filtergRNAs(potential.gRNAs,
                pairOutputFile = pairOutputFile,
                findgRNAsWithREcutOnly = findgRNAsWithREcutOnly,
	        REpatternFile = REpatternFile,
                format = format,  minREpatternSize = minREpatternSize,
                overlap.gRNA.positions = overlap.gRNA.positions)
            REcutDetails  <- gRNAs.RE$gRNAREcutDetails
	    write.table(REcutDetails[order(as.character(
                REcutDetails$ForwardgRNAName)), ], file = REcutDetailFile,
                sep = "\t", row.names = FALSE)
        }
        else if (length(potential.gRNAs) >0)
	{
            gRNAs.RE <- filtergRNAs(potential.gRNAs,
	        findgRNAsWithREcutOnly = findgRNAsWithREcutOnly,
                REpatternFile = REpatternFile, format = format,
                minREpatternSize = minREpatternSize,
                overlap.gRNA.positions = overlap.gRNA.positions)
	    REcutDetails  <- gRNAs.RE$gRNAREcutDetails
	    write.table(REcutDetails[order(as.character(
                REcutDetails$REcutgRNAName)), ], file = REcutDetailFile,
                sep = "\t", row.names = FALSE)
	}
	if (findgRNAsWithREcutOnly)
	{
	    gRNAs  <- gRNAs.RE$gRNAs
        }
	else
	{
	    gRNAs <- potential.gRNAs
	}
        if ( annotatePaired || findPairedgRNAOnly)
	    pairedInformation <- read.table(pairOutputFile, sep = "\t",
                header = TRUE, stringsAsFactors = FALSE)
    }
    else
    {
        if (class(inputFilePath) != "DNAStringSet")
        {
            if (! file.exists(inputFilePath)) {
                stop("inputfile specified as ", inputFilePath, " does not exists!")
            }
            if (format == "fasta" || format == "fastq")
            {
                potential.gRNAs <- readDNAStringSet(inputFilePath, format,
                      use.names = TRUE)
            }
            else
            {
                stop("format needs to be either fasta,fastq or bed!")
            }
        }
        else
        {
            potential.gRNAs <- inputFilePath
            if (length(names(potential.gRNAs)) == 0)
               names(potential.gRNAs) <- paste("gRNAs", 1:length(potential.gRNAs), sep="")
        }
	gRNAs.RE <- filtergRNAs(potential.gRNAs,
            REpatternFile = REpatternFile, format = format,
            minREpatternSize = minREpatternSize,
            overlap.gRNA.positions = overlap.gRNA.positions)
	REcutDetails  <- gRNAs.RE$gRNAREcutDetails
	write.table(
            REcutDetails[order(as.character(REcutDetails$REcutgRNAName)), ],
            file = REcutDetailFile, sep = "\t", row.names = FALSE)
	if (findgRNAsWithREcutOnly)
	{
	    gRNAs  <- gRNAs.RE$gRNAs
	}
	else
	{
	    gRNAs <- potential.gRNAs
	}
	pairedInformation <- ""
    }
    if (length(chromToSearch) == 1 && chromToSearch == "")
    {
	cat("Done. Please check output files in directory ", outputDir, "\n")
        return(gRNAs)
    }
    if (missing(BSgenomeName) || class(BSgenomeName) != "BSgenome") {
        stop("BSgenomeName is required as BSgenome object!")
    }
    if (annotateExon && (missing(txdb) || (class(txdb) != "TxDb" &&
        class(txdb) != "TranscriptDb")))
    {
        stop("To indicate whether an offtarget is inside an exon, txdb is
            required as TxDb object!")
    }
    names(gRNAs) <- gsub( "\t", "", names(gRNAs))
    names(gRNAs) <- gsub( "\n", "", names(gRNAs))
    names(gRNAs) <- gsub( " ", "", names(gRNAs))

    hits <- searchHits2(gRNAs = gRNAs, PAM = PAM, PAM.pattern = PAM.pattern,
        BSgenomeName = BSgenomeName, chromToSearch = chromToSearch,
	chromToExclude = chromToExclude,
        max.mismatch = max.mismatch, PAM.size = PAM.size,
        gRNA.size = gRNA.size, allowed.mismatch.PAM = allowed.mismatch.PAM,
        PAM.location = PAM.location,
        baseEditing = baseEditing, targetBase = targetBase,
        editingWindow = editingWindow.offtargets)
if (dim(hits)[1] > 0)
{
    cat("Building feature vectors for scoring ...\n")
    #save(hits, file = "hits.RData")
    featureVectors <- buildFeatureVectorForScoring(hits = hits,
        canonical.PAM = PAM, gRNA.size = gRNA.size,
        subPAM.position = subPAM.position,
        PAM.location = PAM.location, PAM.size = PAM.size)
    cat("Calculating scores ...\n")
    #save(featureVectors, file="featureVectors.RData")
    if ( scoring.method ==  "CFDscore")
        scores <- getOfftargetScore2(featureVectors,
            subPAM.activity = subPAM.activity,
            mismatch.activity.file = mismatch.activity.file)
    else
        scores <- getOfftargetScore(featureVectors, weights = weights)
    #write.table(scores, file="testScore2.xls", sep="\t", row.names=FALSE)
    cat("Annotating, filtering and generating reports ...\n")
    #saveRDS(scores, file="scores.RDS")
    offTargets <- filterOffTarget(scores = scores, outputDir = outputDir,
        BSgenomeName = BSgenomeName, fetchSequence = fetchSequence, txdb = txdb,
            orgAnn = orgAnn, ignore.strand = ignore.strand,
	    min.score = min.score, topN = topN,
            topN.OfftargetTotalScore = topN.OfftargetTotalScore,
            upstream = upstream, downstream = downstream,
            annotateExon = annotateExon, baseBeforegRNA = baseBeforegRNA,
	    baseAfterPAM = baseAfterPAM, gRNA.size = gRNA.size,
            PAM.location = PAM.location, PAM.size = PAM.size,
            featureWeightMatrixFile = featureWeightMatrixFile,
            rule.set = rule.set, chrom_acc = chrom_acc,
            calculategRNAefficacyForOfftargets = calculategRNAefficacyForOfftargets)
  #saveRDS(offTargets, file = "offTargets.RDS")
    cat("Done annotating\n")
    summary <- read.table(paste(outputDir, "Summary.xls", sep = ""), sep = "\t",
        header = TRUE, stringsAsFactors = FALSE)
    if (dim(summary)[2] == 1)
    	summary <- as.data.frame(t(data.matrix(offTargets$summary)))
    for (i in grep("topOfftarget", names(summary)))
    {
        y <- as.character(summary[,i])
        y[is.na(y)] <- ""
	summary[, i] = y
    }
    if (findgRNAs && (annotatePaired || findPairedgRNAOnly))
    {
        cat("Add paired information...\n")
        PairedgRNAName <- unlist(lapply(1:dim(summary)[1], function(i) {
            as.character(gsub("^\\s+|\\s+$", "",
                paste(unique(pairedInformation[as.character(
                pairedInformation$ForwardgRNAName) == as.character(
                summary$names[i]),]$ReversegRNAName),
                unique(pairedInformation[as.character(
                pairedInformation$ReversegRNAName) == as.character(
                summary$names[i]),]$ForwardgRNAName),
                collapse = " ")))
        }))
    }
    cat("Add RE information...\n")
    if (findPairedgRNAOnly && findgRNAs)
    {
        REname <- unlist(lapply(1:dim(summary)[1], function(i) {
            gsub("^\\s+|\\s+$", "", gsub("NA", "",
                paste(unique(REcutDetails[as.character(
                REcutDetails$ForwardREcutgRNAName) == as.character(
                summary$names[i]),]$ForwardREname),
                unique(REcutDetails[as.character(
                REcutDetails$ReverseREcutgRNAName) ==
                as.character(summary$names[i]), ]$ReverseREname),
                collapse = " ")))
       }))
       summary <- cbind(summary, PairedgRNAName, REname)
    }
    else
    {
        REname <- unlist(lapply(1:dim(summary)[1], function(i) {
            gsub("^\\s+|\\s+$", "", gsub("NA", "", paste(unique(
                REcutDetails[as.character(REcutDetails$REcutgRNAName) ==
                as.character(summary$names[i]), ]$REname), collapse = " ")))
        }))
        summary <- cbind(summary, REname)
    }
	seq <- as.character(summary$gRNAsPlusPAM)
	cat("write gRNAs to bed file...\n")
	on.target <- offTargets$offtargets
	on.target <- unique(subset(on.target,
            on.target$n.mismatch == 0 & on.target$isCanonicalPAM ==1))
	#   as.character(on.target$gRNAPlusPAM) == as.character(on.target$OffTargetSequence)))
	if (dim(on.target)[1] >0)
        {
	   gRNA.bed <- unique(cbind(as.character(on.target$chrom),as.character(on.target$chromStart),
		as.character(on.target$chromEnd), as.character(on.target$name),
		as.numeric(as.character(on.target$gRNAefficacy)) * 1000,
		as.character(on.target$strand),
		as.character(on.target$chromStart),
		as.character(on.target$chromEnd)))
	   if (!useScore)
	   {
		gRNA.bed <- cbind(gRNA.bed, rep("255,0,0",dim(gRNA.bed)[1]))
		gRNA.bed[gRNA.bed[,6] == "-",9] = "0,255,0"
	   }
	#### UCSC genome browser is 0-based instead of 1 based index
	   gRNA.bed[, 2] = as.numeric(gRNA.bed[, 2]) -1
	   gRNA.bed[, 3] = as.numeric(gRNA.bed[, 3])
	   gRNA.bed[gRNA.bed[,6] == "+" ,7] <- as.numeric(gRNA.bed[gRNA.bed[,6] == "+" ,2]) +
		min(overlap.gRNA.positions) - 1
           gRNA.bed[gRNA.bed[,6] == "-" ,7] <- as.numeric(gRNA.bed[gRNA.bed[,6] == "-" ,3]) -
		max(overlap.gRNA.positions)
           gRNA.bed[gRNA.bed[,6] == "+", 8] <- as.numeric(gRNA.bed[gRNA.bed[,6] == "+", 2]) +
		max(overlap.gRNA.positions)
	   gRNA.bed[gRNA.bed[,6] == "-", 8] <- as.numeric(gRNA.bed[gRNA.bed[,6] == "-", 3]) -
		 min(overlap.gRNA.positions) + 1
	   write.table("track name=\"gRNA sites\" description=\"CRISPRseek\" visibility=2 useScore=1 itemRgb=\"On\"", file=bedFile, col.names=FALSE, row.names=FALSE, quote = FALSE)
	  write.table(gRNA.bed, file=bedFile, sep=" ", row.names=FALSE, col.names=FALSE, append=TRUE, quote = FALSE)
	  on.target <- unique(cbind(as.character(on.target$name),
			as.character(on.target$forViewInUCSC),
			as.character(on.target$extendedSequence),
			as.character(on.target$gRNAefficacy)
                        ))
	  colnames(on.target) = c("names", "forViewInUCSC", "extendedSequence", "gRNAefficacy")
	  if (useEfficacyFromInputSeq)
	  {
		on.target <- as.data.frame(on.target[,1:2])
		inputEfficacy <- read.table(efficacyFile, sep="\t", header = TRUE,
			stringsAsFactors=FALSE)
		inputEfficacy <- as.data.frame(cbind(name = inputEfficacy$name,
		        extendedSequence = inputEfficacy$extendedSequence,
			gRNAefficacy = inputEfficacy$gRNAefficacy))
		on.target <- merge(on.target, inputEfficacy, by.x="names", by.y ="name")
	  }
          if(dim(on.target)[1] >0)
	     summary <- unique(merge(on.target, summary, by="names", all=TRUE))
	  write.table(summary[order(as.character(summary$names)), ],
             file = paste(outputDir, "Summary.xls", sep = ""),
             sep = "\t", row.names = FALSE)
	  cat("Scan for REsites in flanking region...\n")
	  if (outputUniqueREs && !missing(BSgenomeName) &&
               class(BSgenomeName) == "BSgenome")
	  {
	    REs.isUnique100 <- uniqueREs(REcutDetails = REcutDetails,
		   summary = summary, offTargets$offtargets, scanUpstream = 100,
		   scanDownstream =100, BSgenomeName = BSgenomeName)
	    REs.isUnique50 <- uniqueREs(REcutDetails = REcutDetails,
		   summary = summary, offTargets$offtargets, scanUpstream = 50,
		   scanDownstream = 50, BSgenomeName = BSgenomeName)
	    summary <- cbind(summary, uniqREin200 = REs.isUnique100,
                uniqREin100 = REs.isUnique50)
            summary$uniqREin200 <- as.character(summary$uniqREin200)
            summary$uniqREin100 <- as.character(summary$uniqREin100)
	  }
	 else
	{
	   REs.isUnique100 = ""
       	   REs.isUnique50 = ""
	}
    }
    else
    {
       warnings("No on-target found for the input gRNAs with your search criteria!")
       gRNA.bed = ""
       REs.isUnique100 = ""
       REs.isUnique50 = ""
    }
    if (foldgRNAs)
     {
        source(system.file("extdata/foldgRNAs.R",package = "CRISPRseek"))
	gRNAs.withoutPAM <- substr(as.character(summary$gRNAsPlusPAM), 1, gRNA.size)
        folded.gRNAs <- foldgRNAs(gRNAs.withoutPAM, gRNA.backbone = gRNA.backbone,
           temperature = temperature)
	if (length(dim(folded.gRNAs)) > 0)
	{
	   if (dim(folded.gRNAs)[1] >1)
	      summary <- cbind(summary, folded.gRNAs[,-1])
	   else
	      summary <- data.frame(c(summary, folded.gRNAs[,-1]))
	}
     }
    #write.table(summary[order(as.character(summary$forViewInUCSC)), ],
    ### even there is no perfect target for a gRNA, it will be kept in the summary file
    ### need to calculate the topN offtarget score and distance correctly yet if include those gRNAs without target

     gRNAs.notInGenome <- setdiff(names(gRNAs), summary$names)
     if (length(gRNAs.notInGenome) > 0)
     {
         dat2 <- data.frame(matrix(nrow = length(gRNAs.notInGenome), ncol = dim(summary)[2]))
         colnames(dat2) <- colnames(summary)
         dat2$names <- gRNAs.notInGenome

         dat2$gRNAsPlusPAM <- paste(substr(as.character(gRNAs[names(gRNAs) %in% gRNAs.notInGenome]), 1, gRNA.size), PAM, sep ="")
         summary <- rbind(summary, dat2)
     }
     if (dim(on.target)[1] == 0)
        write.table(summary[order(as.character(summary$names)), ],
           file = paste(outputDir, "Summary.xls", sep = ""),
           sep = "\t", row.names = FALSE)
     else
        write.table(summary[order(as.character(summary$forViewInUCSC)), ],
           file = paste(outputDir, "Summary.xls", sep = ""),
           sep = "\t", row.names = FALSE)
    if (predIndelFreq) {
        if (predictIndelFreq.onTargetOnly)
		targets <- unique(subset(offTargets$offtargets,
                     offTargets$offtargets$n.mismatch == 0 & offTargets$offtargets$isCanonicalPAM ==1))
        else
		targets <- subset(offTargets$offtargets, offTargets$offtargets$isCanonicalPAM == 1)

        extendedSequence <- getExtendedSequence(targets,
                 BSgenomeName = BSgenomeName,
                 baseBeforegRNA =  baseBeforegRNA.indelFreq,
                 baseAfterPAM = baseAfterPAM.indelFreq, forMethod = method.indelFreq)
        tryCatch((
	    indelFreqFS <- predictRelativeFreqIndels(extendedSequence, method = method.indelFreq)
         ), error = function(e) {print(e); })
        if (exists("indelFreqFS"))
        {
            fs <- unlist(lapply(indelFreqFS, function(x) { x$fs }))
            indelFreq <- lapply(indelFreqFS, function(x) {x$indel})

            entropy <- unlist(lapply(indelFreq, function(x) {
               if (length(x) > 1)
                   sum(-as.numeric(x[,2])/100 *  log(as.numeric(x[,2])/100, base = 450), na.rm = TRUE)
               else
	           NA
            }))
            fs2 <- data.frame(cbind(names = as.character(targets[,1]), frameshift = fs,
               entropy = entropy, n.mismatch = as.character(targets$n.mismatch)))
            fs2[,1] <- as.character(fs2[,1])
            summary <- data.frame(summary)
            summary[,1] <- as.character(summary[,1])

            summary <- merge(subset(fs2, fs2[,4] == 0)[,-4], summary, all.y = TRUE)

            write.table(summary[order(as.character(summary$forViewInUCSC)), ],
               file = paste(outputDir, "Summary.xls", sep = ""),
               sep = "\t", row.names = FALSE)

            names(indelFreq) <- paste(targets[,1], targets[,2], targets[,3],
		 sep= ",")

            if (!predictIndelFreq.onTargetOnly)
            {
                offTargets$offtargets[,3] <- as.character(offTargets$offtargets[,3])
                fs3 <- cbind(OffTargetSequence =  as.character(targets[,3]), frameshift = fs, entropy = entropy)
                targets <- merge(offTargets$offtargets, fs3, all.x = TRUE)
                offTargets$offtargets <- targets
                write.table(targets,  file = paste(outputDir, "OfftargetAnalysis.xls", sep = ""),
                    sep = "\t", row.names = FALSE)
             }
            cat("Done. Please check output files in directory \n", outputDir, "\n")
            list(on.target=on.target, summary=summary, offtarget = offTargets$offtargets,
                 gRNAs.bedFormat=gRNA.bed, REcutDetails = REcutDetails,
                 REs.isUnique100 = REs.isUnique100, REs.isUnique50 = REs.isUnique50,
                 indelFreq = indelFreq, fs2=fs2)
        }
        else
        {
            cat("Done. Please check output files in directory \n", outputDir, "\n")
            list(on.target=on.target, summary=summary, offtarget = offTargets$offtargets,
                 gRNAs.bedFormat=gRNA.bed, REcutDetails = REcutDetails,
                 REs.isUnique100 = REs.isUnique100, REs.isUnique50 = REs.isUnique50)
       }
    }
    else {
        cat("Done. Please check output files in directory \n", outputDir, "\n")
        list(on.target=on.target, summary=summary, offtarget = offTargets$offtargets,
		 gRNAs.bedFormat=gRNA.bed, REcutDetails = REcutDetails,
		 REs.isUnique100 = REs.isUnique100, REs.isUnique50 = REs.isUnique50)
    }
}
else
{
  if (PAM.location == "3prime")
      x <- paste(substr(as.character(gRNAs), 1, gRNA.size), PAM, sep ="")
  else
      x <- paste(PAM, substr(as.character(gRNAs), 1, gRNA.size), sep ="")
  summary <- cbind(names = names(gRNAs), gRNAsPlusPAM = x,top5OfftargetTotalScore = rep("NA", length(gRNAs)),
   	top10OfftargetTotalScore =  rep("NA", length(gRNAs)),
	top1Hit.onTarget.MMdistance2PAM =  rep("NA", length(gRNAs))
      )
  write.table(summary,  file = paste(outputDir, "Summary.xls", sep = ""),
        sep = "\t", row.names = FALSE)
  summary
}
}
