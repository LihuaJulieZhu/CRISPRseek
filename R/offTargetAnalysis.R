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
#' @param inputFilePath Path to an input sequence file or a `DNAStringSet` 
#' object containing sequences to be searched for potential gRNAs.
#' @param format Defaults to "fasta". Format of the input file, "fasta", 
#' "fastq", and "bed" are supported.
#' @param header Defaults to FALSE. Indicates whether the input file contains 
#' header. Only relevant when `format` is set to "bed".
#' @param gRNAoutputName Defaults to "test". Specifies the name of the gRNA 
#' output file when `inputFilePath` is a `DNAStringSet` object instead of 
#' a file path.
#' @param findgRNAs Defaults to TRUE. Specifies whether to find gRNAs from the 
#' sequences in `inputFilePath`. Set to FALSE if the input file already contains 
#' user-selected gRNAs plus PAM.
#' @param exportAllgRNAs Defaults to "both". Indicates whether to output all 
#' potential gRNAs to a file in fasta format, genbank format, or both.
#' @param findgRNAsWithREcutOnly Defaults to TRUE. Specifies whether to search 
#' for gRNAs that overlap with restriction enzyme recognition sites only.
#' @param REpatternFile Path to a file containing restriction enzyme cut 
#' patterns.
#' @param minREpatternSize Defaults to 4. Minimum restriction enzyme recognition
#' pattern length required for the enzyme pattern to be searched for.
#' @param overlap.gRNA.positions Defaults to `c(17, 18)`. Specifies the required 
#' overlapping positions of the gRNA and restriction enzyme cut site. For Cpf1, 
#' you can set it to `c(19, 23)`.
#' @param findPairedgRNAOnly Defaults to FALSE. Specifies whether to search only
#' for paired gRNAs in such an orientation that the first one is on the minus 
#' strand (reverse gRNA) and the second one is on plus strand (forward gRNA).
#' @param annotatePaired Defaults to TRUE. Specifies whether to output paired 
#' gRNA information.
#' @param paired.orientation The "PAMin" orientation refers to the scenario 
#' where the two adjacent PAMs on the sense and antisense strands face inward 
#' toward each other, such as in "N21GG" and "CCN21". In contrast, the "PAMout" 
#' orientation occurs when the PAMs face away from each other, as seen in 
#' "CCN21" and "N21GG".
#' @param min.gap Defaults to 0. Minimum distance between two oppositely 
#' oriented gRNAs to be considered as valid paired gRNAs. 
#' @param enable.multicore Defaults to FALSE. Indicates whether to enable 
#' parallel. For super long sequences with lots of gRNAs, set it to TRUE.
#' @param n.cores.max Defaults to 6. Specifies the maximum number of cores to 
#' use in multicore mode. Set it to 1 to disable multicore processing for 
#' small dataset.
#' @param max.gap Defaults to 20. Specifies the maximum distance between two 
#' oppositely oriented gRNAs to be considered as valid paired gRNAs.
#' @param gRNA.name.prefix Defaults to "gRNA". Specifies the prefix used when 
#' assigning names to detected gRNAs.
#' @param gRNA.size Defaults to 20. The size of the gRNA.
#' @param PAM Defaults to "NGG". Defines the protospacer adjacent motif sequence.
#' @param PAM.size Defaults to `width(PAM)`. Specifies the PAM length.
#' @param PAM.location Defaults to "3prime" (for spCas9). Specifies the PAM 
#' location relative to the protospacer sequence. Set to "5prime" for cpf1
#' because its PAM is located at the 5 prime end of the protospacer.
#' @param PAM.pattern Defaults to "NNG$|NGN$" (for spCas9). Specifies the 
#' regular expression of PAM. For cpf1, set to "^TTTN" since its PAM is at the 5
#' prime end.
#' @param allowed.mismatch.PAM Defaults to 1. Maximum number of mismatches 
#' allowed in the PAM sequence for off-target search. The default value 1 allows
#' "NGN" and "NNG" PAM patterns for off-target identification.
#' @param BSgenomeName A `BSgenome` object containing the target genome 
#' sequence, used for off-target search. Please refer to available genomes in 
#' the "BSgenome" package. For example,
#' \itemize{
#' \item{BSgenome.Hsapiens.UCSC.hg19} - for hg19,
#' \item{BSgenome.Mmusculus.UCSC.mm10} - for mm10
#' \item{BSgenome.Celegans.UCSC.ce6} - for ce6
#' \item{BSgenome.Rnorvegicus.UCSC.rn5} - for rn5
#' \item{BSgenome.Drerio.UCSC.danRer7} - for Zv9
#' \item{BSgenome.Dmelanogaster.UCSC.dm3} - for dm3
#' }
#' @param genomeSeqFile Alternative to `BSgenomeName`. Specifies the path to a 
#' custom target genome file in FASTA format, used for off-target search. It is 
#' applicable when `BSgenomeName` is NOT set. When `genomeSeqFile` is set, the 
#' `annotateExon`, `txdb`, and `orgAnn` parameters will be ignored.
#' @param chromToSearch Defaults to "all", meaning all chromosomes in the target 
#' genome are searched for off-targets. Set to a specific chromosome 
#' (e.g., "chrX") to restrict the search to that chromsome only.
#' @param chromToExclude If set to "", means to search off-targets in 
#' chromosomes specified in `chromToSearch`. By default, to exclude haplotype 
#' blocks from off-target search assuming using `hg19` genome, i.e., 
#' `chromToExclude = c("chr17_ctg5_hap1", "chr4_ctg9_hap1", "chr6_apd_hap1", 
#' "chr6_cox_hap2", "chr6_dbb_hap3", "chr6_mann_hap4", "chr6_mcf_hap5",
#' "chr6_qbl_hap6", "chr6_ssto_hap7")`.
#' @param max.mismatch Defaults to 3. Maximum number of mismatches allowed in 
#' off-target search. Warning: search will be considerably slower if set to a 
#' value greater than 3.
#' @param findOffTargetsWithBulge Defaults to FALSE. Specifies whether to search
#' for off-targets with bulges.
#' @param method.findOffTargetsWithBulge Only applicable if 
#' `findOffTargetsWithBulge = TRUE`. Choose from `c("CasOFFinder_v3.0.0b3")`.
#' @param DNA_bulge Defaults to 2. Maximum number of DNA bulges allowed in 
#' off-target search.
#' @param RNA_bulge Defaults to 2. Maximum number of RNA bulges allowed in 
#' off-target search.
#' @param gRNA.pattern Defaults to NULL (meaning no restriction). Specifies 
#' regular expression or IUPAC Extended Genetic Alphabet to represent gRNA 
#' pattern. E.g. to specify that the gRNA must start with "GG", set it to "^GG".
#' Type `?translatePattern` for a list of IUPAC Extended Genetic Alphabet.
#' @param baseEditing Defaults to FALSE. Specifies whether to design gRNAs for 
#' base editing. If set to TRUE, please set `targetBase` and `editingWidow`.
#' @param targetBase Defaults to "C" (for converting C to T in the CBE system). 
#' Applicable only when `baseEditing = TRUE`. Specifies the target base for base
#' editing systems. Please change it to "A" if you intend to use the ABE system.
#' @param editingWindow Defaults to `4:8` (for the CBE system). Applicable only 
#' when `baseEditing = TRUE`, and specifies the effective editing window. Please 
#' change it accordingly if the system you use have a different editing window.
#' @param editingWindow.offtargets Defaults to `4:8` (for the original CBE
#' system, 1 means the most distal site from the 3' PAM, the most proximal site 
#' from the 5' PAM). Applicable only when `baseEditing = TRUE`. Indicates the 
#' effective editing window to consider for the off-targets search only. Please 
#' change it accordingly if the system you use have a different editing window, 
#' or if you would like to include off-targets with the target base in a larger 
#' editing window.
#' @param primeEditing Defaults to FALSE. Specifies whether to design gRNAs for 
#' prime editing. If set to TRUE, please set `PBS.length`, `RT.template.length`,
#' `RT.template.pattern`, `targeted.seq.length.change`, `bp.after.target.end`,
#' `target.start`, `target.end`, and `corrected.seq` accordingly.
#' @param PBS.length Applicable only when `primeEditing = TRUE`. Specifies the 
#' number of bases to output for primer binding site.
#' @param RT.template.length Defaults to `8:18`. Applicable only when 
#' `primeEditing = TRUE`. Specifies the number of bases required for RT 
#' template. Increase the length if the edit involves a large insertion. Only 
#' gRNAs with a calculated `RT.template.length` within the specified range will 
#' be included in the output. It is calculated as the following:
#' `RT.template.length = target.start – cut.start + (target.end - target.start) +
#' targeted.seq.length.change + bp.after.target.end`.
#' @param RT.template.pattern Defaults to not end with C (per 
#' https://doi.org/10.1038/s41586-019-1711-4). Applicable only when
#' `primeEditing = TRUE`. Specifies the RT template sequence pattern.
#' @param corrected.seq Applicable only when `primeEditing = TRUE`. Specifies 
#' the mutated or inserted sequences after successful editing.
#' @param targeted.seq.length.change Applicable only when `primeEditing = TRUE`.
#' Specifies the change in the targeted sequence length. Set it to 0 for base 
#' changes, positive numbers for insertions, and negative number for deletions. 
#' For example, 10 indicates that the corrected sequence will have a 10-bp 
#' insertion, -10 means that the corrected sequence will have a 10-bp deletion, 
#' and 0 means that only base changes with no change in sequence length.
#' @param bp.after.target.end Defaults to 15. Applicable only when 
#' `primeEditing = TRUE`. Specifies the number of bases to add after the target 
#' change end site as part of the RT template. Refer to `RT.template.length` for
#' how this parameter affects the calculation of `RT.template.length`, which is 
#' used as a filtering criterion during pregRNA selection.
#' @param target.start Defaults to 20. Applicable only when 
#' `primeEditing = TRUE`. Specifies the start location in the input sequence to 
#' make changes, which will be used to obtain the RT template sequence. Refer to
#' `RT.template.length` for how this parameter affects the `RT.template.length`
#' calculation, which is used as a filtering criteria in pregRNA selection.
#' @param target.end Defaults to 20. Applicable only when `primeEditing = TRUE`.
#' Specifies the end location in the input sequence to make changes, which will 
#' be used to obtain the RT template sequence. Refer to `RT.template.length` for
#' how this parameter affects the `RT.template.length` calculation, which is
#' used as a filtering criteria in pregRNA selection.
#' @param primeEditingPaired.output Defaults to "pairedgRNAsForPE.xls". 
#' Applicable only when `primeEditing = TRUE`. Specifies the file path where the
#' pegRNA, second gRNA wit PBS, RT.template, and gRNA sequences will be saved.
#' @param min.score Defaults to 0. Specifies the minimum score of an off-target 
#' to be included in the final output.
#' @param topN Defaults to 1000. Specifies the top N off-targets to be included 
#' in the final output
#' @param topN.OfftargetTotalScore Defaults to 10. Specifies the top N 
#' off-targets used to calculate the total off-target score.
#' @param annotateExon Defaults to TRUE. Specifies whether to indicate if the 
#' off-target is located within an exon.
#' @param txdb A `TxDb` object containing organism-specific annotation data, 
#' required for `annotateExon`. For creating and using a `TxDb` object, refer to
#' the `GenomicFeatures` package. For a list of existing `TxDb` objects, search
#' for annotation packages starting with "Txdb" at
#' http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData,
#' such as
#' \itemize{
#' \item{TxDb.Rnorvegicus.UCSC.rn5.refGene} - for rat
#' \item{TxDb.Mmusculus.UCSC.mm10.knownGene} - for mouse
#' \item{TxDb.Hsapiens.UCSC.hg19.knownGene} - for human
#' \item{TxDb.Dmelanogaster.UCSC.dm3.ensGene} - for Drosophila
#' \item{TxDb.Celegans.UCSC.ce6.ensGene} - for C.elegans
#' }
#' @param orgAnn An `OrgDb` object containing organism-specific annotation 
#' mapping information, required for `annotateExon`.
#' @param ignore.strand Defaults to TRUE. Specifies if strandness should be 
#' ignored when annotating off-targets to genes.
#' @param outputDir Defaults to the current working directory. Specifies the 
#' path to the directory where the analysis results will be saved.
#' @param fetchSequence Defaults to TRUE. Specifies whether to fetch flanking
#' sequences for off-targets.
#' @param upstream Defaults to 200. Specifies the upstream offset from the 
#' off-target start.
#' @param downstream Defaults to 200. Specifies the downstream offset from the 
#' off-target end.
#' @param weights Defaults to `c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 
#' 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 
#' 0.583)` (used in Hsu et al., 2013 cited in the reference section). Applicable
#' only when `scoring.method = Hus-Zhang`. Specifies a numeric vector with a 
#' length equal to the size of the gRNA, containing the corresponding weight 
#' values.
#' @param baseBeforegRNA Defaults to 4. Specifies the number of bases preceding 
#' the gRNA. It is used to calculate gRNA efficiency. Note that for PAMs located 
#' at the 5 prime end, the number of bases should include both the bases before 
#' the PAM sequence and the PAM size.
#' @param baseAfterPAM Defaults to 3 (for spCas9). Specifies the number of bases
#' after PAM. It is used to calculate gRNA efficiency. Note that for PAMs 
#' located on the 5 prime end, the number should include the length of the gRNA
#' plus the extended sequence on the 3 prime end.
#' @param featureWeightMatrixFile By default, the DoenchNBT2014 weight matrix is 
#' used. Specifies the feature weight matrix file used for calculating gRNA 
#' efficiency. To use an alternative matrix, provide a CSV where the first 
#' column contains the significant features and the second column contains the
#' corresponding weights. For details, refer to Doench et al., 2014.
#' @param useScore Defaults to TRUE. Displays in grayscale, with darkness
#' indicating gRNA efficacy. The taller bar represents the Cas9 cutting site.
#' If set to False, efficacy will not be shown. Instead, gRNAs on the plus 
#' strand will be colored red, and gRNAs on the minus strand will be colored 
#' green.
#' @param useEfficacyFromInputSeq Defaults to FALSE. If TRUE, the summary file
#' will contain gRNA efficacy calculated from the input sequences instead of 
#' from off-target analysis. Set it to TRUE if the input sequence is from a 
#' species different from the one used for off-target analysis.
#' @param outputUniqueREs Defaults to TRUE. If set to TRUE, summary file will
#' contain REs unique to the cleavage site within 100 or 200 bases surrounding
#' the gRNA sequence.
#' @param foldgRNAs Defaults to FALSE. If set to TRUE, summary file will contain
#' minimum free energy of the secondary structure of gRNA with gRNA backbone
#' from `GeneRfold` package given that `GeneRfold` package has been installed.
#' @param gRNA.backbone Defaults to the sequence in Sp gRNA backbone. Applicable
#' only when `foldgRNAs = TRUE`. Specifies the gRNA backbone constant region 
#' sequence.
#' @param temperature Defaults to 30. Applicable only when `foldgRNAs = TRUE`.
#' Specifies the temperature in Celsius.
#' @param overwrite Defaults to FALSE. Specifies whether to overwrite the 
#' existing files in the output directory.
#' @param scoring.method Defaults to "Hsu-Zhang". Specifies the method to use 
#' for off-target cleavage rate estimation. Choose from "Hsu-Zhang" and 
#' "CFDscore"
#' @param subPAM.activity Defaults to "hash(AA = 0, AC = 0, AG = 0.259259259, 
#' AT = 0, CA = 0, CC = 0, CG = 0.107142857, CT = 0, GA = 0.069444444, 
#' GC = 0.022222222, GG = 1, GT = 0.016129032, TA = 0, TC = 0, TG = 0.038961039,
#' TT = 0)". Applicable only when `scoring.method = CFDscore`. Specifies a hash 
#' that represents the cleavage rate for each alternative sub PAM sequence 
#' relative to preferred PAM sequence.
#' @param subPAM.position Defaults to `c(22, 23)` (For spCas9 with 20-bp gRNA and 
#' NGG as preferred PAM). Applicable only when `scoring.method = CFDscore`. 
#' Specifies the start and end positions of the sub PAM. For Cpf1, it should be
#' `c(1,2)`.
#' @param rule.set Defaults to "Root_RuleSet1_2014". Specifies a rule set 
#' scoring system for calculating gRNA efficacy. Note that "Root_RuleSet2_2016"
#' requires the following packages with specified version: python 2.7,
#' scikit-learn 0.16.1, pickle, pandas, numpy, and scipy.
#' @param chrom_acc Specifies an optional binary variable indicating chromatin 
#' accessibility information with 1 representing accessible and 0 representing 
#' inaccessible.
#' @param calculategRNAefficacyForOfftargets Defaults to TRUE. Specifies whether
#' to output gRNA efficacy for both off-targets and on-targets. Set to FALSE if
#' only on-target gRNA efficacy is needed to speed up the analysis. For 
#' potential use cases of off-target efficacies, refer to 
#' https://support.bioconductor.org/p/133538/#133661.
#' @param mismatch.activity.file Defaults to use the supplementary Table 19 from
#' Doench et al., Nature Biotechnology 2016. Applicable only when 
#' `scoring.method = CFDscore`. Specifies a CSV file containing the cleavage 
#' rates for all possible types of single nucleotide mismatches at each position
#' of the gRNA.
#' @param predIndelFreq Defaults to FALSE. Specifies whether to output the 
#' predicted INDELs and their frequencies.
#' @param predictIndelFreq.onTargetOnly Defaults to TRUE. Specifies whether to
#' predict INDELs and their frequencies for on-targets only. Typically,
#' researchers are only interested in predicting editing outcome for on-targets,
#' as editing in off-targets is undesirable. Set to FALSE if you want to predict 
#' INDELs and their frequencies for off-targets as well. Note that this will 
#' increase the run time.
#' @param method.indelFreq Defaults to "Lindel". Applicable only when 
#' `predIndelFreq = TRUE`. Specifies the method to be used for predicting 
#' INDELs. Currently, only "Lindel" is supported, though additional methods can 
#' be added upon request. Type `?predictRelativeFreqIndels` to learn more.
#' @param baseBeforegRNA.indelFreq Defaults to 13. Applicable only when
#' `predIndelFreq = TRUE`.
#' @param baseAfterPAM.indelFreq Defaults to 24. Applicable only when
#' `predIndelFreq = TRUE`.
#' @return Four Excel files are generated in the output directory:
#' \item{Summary.xlsx}{ - Summary of the gRNAs}
#' \item{OfftargetAnalysis.xlsx}{ - Detailed information on off-targets}
#' \item{REcutDetails.xlsx}{ - Restriction enzyme cut sites for each gRNA}
#' \item{pairedgRNAs.xlsx}{ - Potential paired gRNAs}
#' @note %% ~~further notes~~
#' @author Lihua Julie Zhu, Kai Hu
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
#' # Load required libraries
#' library(CRISPRseek)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(org.Hs.eg.db)
#' 
#' # Example 1: given FASTA input, search gRNAs and off-targets
#' outputDir <- tempdir()
#' inputFilePath <- system.file("extdata/inputseq.fa", package = "CRISPRseek")
#' REpatternFile <- system.file("extdata/NEBenzymes.fa", package = "CRISPRseek")
#' 
#' results <- offTargetAnalysis(inputFilePath, 
#'                              findPairedgRNAOnly = FALSE,
#'                              findgRNAsWithREcutOnly = TRUE,
#'                              REpatternFile = REpatternFile, 
#'                              annotatePaired = FALSE,
#'                              BSgenomeName = Hsapiens, 
#'                              chromToSearch = "chrX",
#'                              txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                              orgAnn = org.Hs.egSYMBOL, 
#'                              max.mismatch = 1,
#'                              outputDir = outputDir, 
#'                              overwrite = TRUE)
#' 
#' # Example 2: also predict indels and frequecies at target sites
#' results <- offTargetAnalysis(inputFilePath,
#'                              predIndelFreq = TRUE, 
#'                              predictIndelFreq.onTargetOnly = TRUE,
#'                              findgRNAsWithREcutOnly = TRUE,
#'                              findPairedgRNAOnly = FALSE,
#'                              annotatePaired = FALSE,
#'                              BSgenomeName = Hsapiens, 
#'                              chromToSearch = "chrX",
#'                              txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                              orgAnn = org.Hs.egSYMBOL, 
#'                              max.mismatch = 1,
#'                              outputDir = outputDir, 
#'                              overwrite = TRUE)
#' names(results$indelFreq)
#' head(results$indelFreq[[1]])
#'   # Save the indel frequences to tab delimited files, 
#'   # one file for each target or offtarget site.
#' mapply(openxlsx::write.xlsx, results$indelFreq, 
#'        file = paste0(names(results$indelFreq), '.xlsx'))
#' 
#' # Example 3: predict gRNA efficacy using CRISPRscan
#' featureWeightMatrixFile <- system.file("extdata/Morenos-Mateo.csv", 
#'                                        package = "CRISPRseek")
#' 
#' results <- offTargetAnalysis(inputFilePath, 
#'                              rule.set = "CRISPRscan",
#'                              findgRNAsWithREcutOnly = TRUE,
#'                              REpatternFile = REpatternFile, 
#'                              findPairedgRNAOnly = FALSE,
#'                              annotatePaired = FALSE,
#'                              BSgenomeName = Hsapiens, 
#'                              chromToSearch = "chrX",
#'                              txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                              orgAnn = org.Hs.egSYMBOL, 
#'                              max.mismatch = 1,
#'                              baseBeforegRNA = 6, 
#'                              baseAfterPAM = 6,
#'                              featureWeightMatrixFile = featureWeightMatrixFile,
#'                              outputDir = outputDir, 
#'                              overwrite = TRUE)
#' 
#' # Example 4: when PAM is on the 5 prime side, e.g., Cpf1
#' if (interactive()) {
#'   results <- offTargetAnalysis(inputFilePath = 
#'                                system.file("extdata/cpf1-2.fa", 
#'                                            package = "CRISPRseek"), 
#'                                PAM.location = "5prime",
#'                                rule.set = "DeepCpf1",
#'                                PAM.size = 4,
#'                                PAM = "TTTN", 
#'                                PAM.pattern = "^TNNN", 
#'                                findgRNAsWithREcutOnly =  FALSE,
#'                                findPairedgRNAOnly = FALSE,
#'                                annotatePaired = FALSE,
#'                                BSgenomeName = Hsapiens,
#'                                chromToSearch = "chr8",
#'                                txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                                orgAnn = org.Hs.egSYMBOL, max.mismatch = 4,
#'                                baseBeforegRNA = 8, baseAfterPAM = 26,
#'                                overlap.gRNA.positions = c(19, 23),
#'                                useEfficacyFromInputSeq = FALSE,
#'                                outputDir = outputDir,
#'                                overwrite = TRUE, 
#'                                allowed.mismatch.PAM = 2,
#'                                subPAM.position = c(1, 2))
#' }
#' 
#' # Example 5: when PAM is on the 5 prime side, and using Root_RuleSet1_2014
#' results <- offTargetAnalysis(inputFilePath, 
#'                              PAM.location = "5prime",
#'                              PAM = "TGT", 
#'                              PAM.pattern = "^T[A|G]N", 
#'                              findgRNAsWithREcutOnly =  FALSE,
#'                              REpatternFile = REpatternFile, 
#'                              findPairedgRNAOnly = FALSE,
#'                              annotatePaired = FALSE,
#'                              BSgenomeName = Hsapiens, 
#'                              chromToSearch = "chrX",
#'                              txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                              orgAnn = org.Hs.egSYMBOL, 
#'                              max.mismatch = 4,
#'                              outputDir = outputDir, 
#'                              overwrite = TRUE, 
#'                              allowed.mismatch.PAM = 2, 
#'                              subPAM.position = c(1, 2), 
#'                              baseEditing = TRUE, 
#'                              editingWindow = 20, 
#'                              targetBase = "G")
#' 
#' # Example 6: base editor
#' results <- offTargetAnalysis(inputFilePath, 
#'                              baseEditing = TRUE,
#'                              editingWindow = 10:20, 
#'                              targetBase = "A",
#'                              findgRNAsWithREcutOnly = FALSE,
#'                              REpatternFile = REpatternFile, 
#'                              findPairedgRNAOnly = FALSE,
#'                              annotatePaired = FALSE,
#'                              BSgenomeName = Hsapiens, 
#'                              chromToSearch = "chrX",
#'                              txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                              orgAnn = org.Hs.egSYMBOL, 
#'                              max.mismatch = 4,
#'                              PAM.location = "5prime",
#'                              PAM = "TGT", 
#'                              PAM.pattern = "^T[A|G]N", 
#'                              allowed.mismatch.PAM = 2,
#'                              subPAM.position = c(1, 2),
#'                              outputDir = outputDir, 
#'                              overwrite = TRUE)
#' 
#' # Example 7: prime editor
#' inputFilePath <- DNAStringSet(paste0("CCAGTTTGTGGATCCTGCTCTGGTGTCCTCCACACC",
#'                                      "AGAATCAGGGATCGAAAACTCATCAGTCGATGCGAG", 
#'                                      "TCATCTAAATTCCGATCAATTTCACACTTTAAACG"))
#' results <- offTargetAnalysis(inputFilePath,
#'                              primeEditing = TRUE, 
#'                              overlap.gRNA.positions = c(17, 18),
#'                              PBS.length = 15,
#'                              corrected.seq = "T",
#'                              RT.template.pattern = "D$",
#'                              RT.template.length = 8:30,
#'                              targeted.seq.length.change = 0,
#'                              bp.after.target.end = 15,
#'                              target.start = 20,
#'                              target.end = 20,
#'                              paired.orientation = "PAMin",
#'                              findPairedgRNAOnly = TRUE,
#'                              BSgenomeName = Hsapiens, 
#'                              chromToSearch = "chrX",
#'                              txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                              orgAnn = org.Hs.egSYMBOL, 
#'                              max.mismatch = 1,
#'                              outputDir = outputDir, 
#'                              overwrite = TRUE,
#'                              PAM.size = 3,
#'                              gRNA.size = 20,
#'                              min.gap = 20, 
#'                              max.gap = 90)
#'
#' @importFrom hash hash
#' @importFrom utils read.csv write.table read.table
#' @importFrom GenomicRanges intersect setdiff
#' @importFrom Biostrings writeXStringSet readDNAStringSet
#' @importFrom BiocGenerics rbind as.data.frame cbind unlist lapply
#' @importFrom rlang abort inform arg_match call_match warn
#' @export
offTargetAnalysis <- function(inputFilePath = NULL, 
                              format = c("fasta", "fastq", "bed"), 
                              header = FALSE,
                              gRNAoutputName = "test", 
                              findgRNAs = TRUE,
                              exportAllgRNAs = c("all", "fasta", "genbank", "no"),
                              findgRNAsWithREcutOnly = FALSE,
                              REpatternFile = REpatternFile_default(),
                              minREpatternSize = 4,
                              overlap.gRNA.positions = c(17, 18), 
                              findPairedgRNAOnly = FALSE,
                              annotatePaired = TRUE, 
                              paired.orientation = c("PAMout", "PAMin"),
                              enable.multicore = FALSE, 
                              n.cores.max = 6,
                              min.gap = 0, 
                              max.gap = 20, 
                              gRNA.name.prefix = NULL,
                              gRNA.size = 20, 
                              PAM = "NGG", 
                              PAM.size = width(PAM),
                              PAM.pattern = "NNG$|NGN$",
                              BSgenomeName = NULL,
                              genomeSeqFile = NULL,
                              chromToSearch = "all",
                              chromToExclude = chromToExclude_default,
                              max.mismatch = 3,
                              allowed.mismatch.PAM = 1,
                              gRNA.pattern = NULL,
                              baseEditing = FALSE, 
                              targetBase = "C", 
                              editingWindow = 4:8,
                              editingWindow.offtargets = 4:8,
                              primeEditing = FALSE,
                              PBS.length = 13L,
                              RT.template.length = 8:28,
                              RT.template.pattern = "D$",
                              corrected.seq = NULL,
                              targeted.seq.length.change = NULL,
                              bp.after.target.end = 15L,
                              target.start = NULL,
                              target.end = NULL,
                              primeEditingPaired.output = "pairedgRNAsForPE.xls",
                              min.score = 0, 
                              topN = 1000,
                              topN.OfftargetTotalScore = 10,
                              annotateExon = TRUE, 
                              txdb = NULL, 
                              orgAnn = NULL, 
                              ignore.strand = TRUE, 
                              outputDir = getwd(),
                              fetchSequence = TRUE, 
                              upstream = 200, 
                              downstream = 200,
                              weights = weights_default,
                              baseBeforegRNA = 4, 
                              baseAfterPAM = 3,
                              featureWeightMatrixFile = featureWeightMatrixFile_default(),
                              useScore = TRUE, 
                              useEfficacyFromInputSeq = FALSE,
                              outputUniqueREs = TRUE, 
                              foldgRNAs = FALSE,
                              gRNA.backbone = gRNA.backbone_default,
                              temperature = 37,
                              overwrite = FALSE,
                              scoring.method = c("Hsu-Zhang", "CFDscore"),
                              subPAM.activity = subPAM.activity_default,
                              subPAM.position = c(22, 23),
                              PAM.location = "3prime",
                              rule.set = c("Root_RuleSet1_2014", "Root_RuleSet2_2016", "CRISPRscan", "DeepCpf1"),
                              chrom_acc = NULL,
                              calculategRNAefficacyForOfftargets = TRUE,
                              mismatch.activity.file = mismatch.activity.file_default(),
                              predIndelFreq = FALSE,
                              predictIndelFreq.onTargetOnly = TRUE,
                              method.indelFreq = "Lindel",
                              baseBeforegRNA.indelFreq = 13,
                              baseAfterPAM.indelFreq = 24,
                              findOffTargetsWithBulge = FALSE,
                              method.findOffTargetsWithBulge = c("CasOFFinder_v3.0.0b3"),
                              DNA_bulge = 2,
                              RNA_bulge = 2) {
  
  ### Step0: params validation: start ###
  inform("Validating arguments ...")
  # Check if BSgenomeName or genomeSeqFile is specified:
  if (is.null(BSgenomeName) && is.null(genomeSeqFile)) {
    if (!chromToSearch == "") {
      abort("Please specify either BSgenomeName (BSgenome format) or genomeSeqFile (FASTA format) to search for off-targets.")
    }
  } else if (!is.null(BSgenomeName) && !is.null(genomeSeqFile)) {
    abort("Both BSgenomeName and genomeSeqFile are provided, only one is allowed!")
  } else if (!is.null(genomeSeqFile)) {
    if (!file.exists(genomeSeqFile)) {
      abort(paste0("File ", genomeSeqFile, " does not exist!"))
    }
    if (annotateExon || !is.null(txdb) || !is.null(BSgenomeName) || fetchSequence) {
      warn("Since genomeSeqFile is supplied, annotateExon, BSgenomeName, txdb, and featchSequence will be ignored!")
      annotateExon <- FALSE
      txdb <- NULL
      BSgenomeName <- NULL
      fetchSequence <- FALSE
    }
  } else if (!inherits(BSgenomeName, "BSgenome")) {
    abort("BSgenomeName must be a BSgenome object!")
  }

  # Check if inputFilePath is is.null:
  if (is.null(inputFilePath)) {
    abort("inputFilePath containing the searching sequence, coordinate or a DNAStringSet object is required!")
  }
  
  # Match arguments:
  format <- arg_match(format)
  if (format == "bed") {
    if (is.null(BSgenomeName) || class(BSgenomeName) != "BSgenome") {
      abort("BSgenomeName is required as BSgenome object when input file is in bed format!")
    }
    inputFilePath <- getSeqFromBed(inputFilePath, header = header, BSgenomeName = BSgenomeName)
    #### format for filtergRNAs
    format <- "fasta"
  }
  scoring.method <- arg_match(scoring.method)
  exportAllgRNAs <- arg_match(exportAllgRNAs)
  rule.set <- arg_match(rule.set)
  paired.orientation <- arg_match(paired.orientation)
  if (findOffTargetsWithBulge) {
    method.findOffTargetsWithBulge <- arg_match(method.findOffTargetsWithBulge)
    if (method.findOffTargetsWithBulge == "CasOFFinder_v3.0.0b3") {
      method.findOffTargetsWithBulge = "3.0.0b3"
    }
  }
  
  # Convert empty arguments to NULL:
  if (length(chromToSearch) == 1 && chromToSearch == "") {
    chromToSearch <- NULL
  }
  if (length(chromToExclude) == 1 && chromToExclude == "") {
    chromToSearch <- NULL
  }
  
  # Fix PAM.pattern:
  PAM.pattern <- fixPAMpattern(PAM.pattern, PAM.location)
  
  # Ensure at least DNA.bulge or RNA.bulge > 0 if find offtarget with bulge:
  if (findOffTargetsWithBulge) {
    if (DNA_bulge <= 0 && RNA_bulge <= 0) {
      abort("Either DNA_bulge or RNA_bulge must be greater than 0!")
    }
  }
  
  # Ensure the length of PAM, PAM.pattern, and PAM.size are consistent:
  checkPAMsize(PAM, PAM.size)
  
  # Fix parameters per rule.set:
  tem <- fixRuleSet(rule.set, scoring.method, subPAM.activity, PAM.location, PAM.size, gRNA.size)
  baseBeforegRNA <- tem[[1]]
  baseAfterPAM <- tem[[2]]
  subPAM.activity <- tem[[3]]

  # Validate parameters per scoring.method:
  validateScoringMethod(scoring.method, mismatch.activity.file, weights, gRNA.size)
  
  # Check REpattern file:
  checkREpatternFile(findgRNAsWithREcutOnly, findgRNAs, REpatternFile)
  
  # Check required libraries:
  checkDependency(foldgRNAs = foldgRNAs)
  
  # Prepare outputDir:
  prepOutputDir(outputDir, overwrite)   
  
  # Prepare other files:
  tem <- prepOtherFiles(annotatePaired, findPairedgRNAOnly, inputFilePath, outputDir, gRNAoutputName)
  pairOutputFile <- tem[[1]]
  REcutDetailFile <- tem[[2]]
  bedFile <- tem[[3]]
  gRNAoutputName <- tem[[4]]
  inform("Validating arguments: done!")
  ### Step0: params validation: finish ###
  
  ### Step0b: prepare lists of grouped arguments: start ###
  # efficacyFile <- tempfile()
  efficacyFile <- file.path(outputDir, "gRNAefficacy.xlsx") # must be hard-coded since other functions rely on it
  arg_groups <- getArgGroups(inputFilePath = inputFilePath,
                             overlap.gRNA.positions = overlap.gRNA.positions,
                             baseEditing = baseEditing, 
                             targetBase = targetBase, 
                             editingWindow = editingWindow,
                             editingWindow.offtargets = editingWindow.offtargets,
                             primeEditing = primeEditing,
                             findPairedgRNAOnly = findPairedgRNAOnly,
                             annotatePaired = annotatePaired,
                             paired.orientation = paired.orientation,
                             pairOutputFile = pairOutputFile,
                             PAM = PAM,
                             PAM.location = PAM.location,
                             PAM.size = PAM.size,
                             gRNA.pattern = gRNA.pattern,
                             gRNA.size = gRNA.size,
                             min.gap = min.gap,
                             max.gap = max.gap,
                             name.prefix = gRNA.name.prefix,
                             format = format,
                             rule.set = rule.set,
                             chrom_acc = chrom_acc,
                             PBS.length = PBS.length,
                             RT.template.length = RT.template.length,
                             RT.template.pattern = RT.template.pattern,
                             targeted.seq.length.change = targeted.seq.length.change,
                             bp.after.target.end = bp.after.target.end,
                             target.start = target.start,
                             target.end = target.end,
                             primeEditingPaired.output = primeEditingPaired.output,
                             corrected.seq = corrected.seq,
                             featureWeightMatrixFile = featureWeightMatrixFile,
                             baseBeforegRNA = baseBeforegRNA,
                             baseAfterPAM = baseAfterPAM,  
                             efficacyFile = efficacyFile,
                             enable.multicore = enable.multicore,
                             n.cores.max = n.cores.max,
                             BSgenomeName = BSgenomeName, # searchHit2
                             genomeSeqFile = genomeSeqFile, # searchHit
                             txdb = txdb,
                             chromToSearch = chromToSearch, #
                             chromToExclude = chromToExclude, #
                             max.mismatch = max.mismatch, #
                             PAM.pattern = PAM.pattern, #
                             allowed.mismatch.PAM = allowed.mismatch.PAM, #
                             outputDir = outputDir, ## outputPotentialgRNAs
                             exportAllgRNAs = exportAllgRNAs, ##
                             gRNAoutputName = gRNAoutputName, ##
                             findgRNAs = findgRNAs, # add paired and RE info to summary, getPotentialgRNA and filtergRNA can't use it
                             findgRNAsWithREcutOnly = findgRNAsWithREcutOnly, ### filtergRNAs
                             REpatternFile = REpatternFile, ###
                             minREpatternSize = minREpatternSize, ###
                             subPAM.position = subPAM.position, #### getOfftargetScoreWrap 
                             subPAM.activity = subPAM.activity, ####
                             scoring.method = scoring.method, ####
                             mismatch.activity.file = mismatch.activity.file, ####
                             weights = weights, ####
                             fetchSequence = fetchSequence, ##### getOfftargetSummary
                             orgAnn = orgAnn, #####
                             ignore.strand = ignore.strand, #####
                             min.score = min.score, #####
                             topN = topN, #####
                             topN.OfftargetTotalScore = topN.OfftargetTotalScore, #####
                             upstream = upstream, #####
                             downstream = downstream, #####
                             annotateExon = annotateExon, #####
                             calculategRNAEfficacy = TRUE,
                             calculategRNAefficacyForOfftargets = calculategRNAefficacyForOfftargets, #####
                             DNA_bulge = DNA_bulge, ###### getOfftargetWithBulge
                             RNA_bulge = RNA_bulge, ######
                             method.findOffTargetsWithBulge = method.findOffTargetsWithBulge ######
                             )
  findgRNAs_args_core <- arg_groups[[1]]
  findgRNAs_args_prime <- arg_groups[[2]]
  findgRNAs_args_efficacy <- arg_groups[[3]]
  searchHits_args <- arg_groups[[4]]
  outputPotentialgRNAs_args <- arg_groups[[5]]
  filtergRNAs_args <- arg_groups[[6]]
  getgRNASummary_args <- arg_groups[[7]]
  getOfftargetScoreWrap_args <- arg_groups[[8]] 
  getOfftargetSummary_args <- arg_groups[[9]]
  getOfftargetWithBulge_args <- arg_groups[[10]]
  filterCasOffinder_args <- arg_groups[[11]]
  ### Step0b: prepare lists of grouped arguments: finish ###
  
  ### Step1: find and output potential gRNAs: start ###
  potential.gRNAs <- getPotentialgRNAs(findgRNAs = findgRNAs,
                                       chromToSearch = chromToSearch, 
                                       inputFilePath = inputFilePath,
                                       format = format,
                                       useEfficacyFromInputSeq = useEfficacyFromInputSeq, 
                                       findgRNAs_args_core = findgRNAs_args_core, 
                                       findgRNAs_args_prime = findgRNAs_args_prime, 
                                       findgRNAs_args_efficacy = findgRNAs_args_efficacy,
                                       outputPotentialgRNAs_args = outputPotentialgRNAs_args)
  if (length(potential.gRNAs) == 0) {
    return(message("no gRNAs found!"))
  }
  ### Step1: find and output potential gRNAs: finish ###

  ### Step2: filter and output potential gRNAs: start ###
  tem_res <- getFiltergRNAs(potential.gRNAs = potential.gRNAs,
                            findgRNAs = findgRNAs,
                            filtergRNAs = filtergRNAs, 
                            filtergRNAs_args = filtergRNAs_args, 
                            findgRNAsWithREcutOnly = findgRNAsWithREcutOnly, 
                            annotatePaired = annotatePaired,
                            findPairedgRNAOnly = findPairedgRNAOnly, 
                            pairOutputFile = pairOutputFile,
                            REcutDetailFile = REcutDetailFile)
  gRNAs <- tem_res[[1]]
  pairedInformation <- tem_res[[2]]
  REcutDetails <- tem_res[[3]]
  ### Step2: filter and output potential gRNAs: finish ###
  
  ### Step3: search for off targets: start ###
  if (is.null(chromToSearch)) {
    inform(paste0("All done! Please check output files in the directory: ", outputDir))
    return(gRNAs)
  } else if (annotateExon && (is.null(txdb) || !inherits(txdb, c("TxDb", "TranscriptDb")))) {
    stop("To indicate whether an off-target is inside an exon, txdb is required as a TxDb object!")
  }
  inform("Searching for off-targets ...")
    ## 3-1: Without bulge:
  if (is.null(genomeSeqFile)) {
    hits <- do.call(searchHits2, c(list(gRNAs = gRNAs, BSgenomeName = BSgenomeName), searchHits_args))
  } else {
    hits <- do.call(searchHitsFa, c(list(gRNAs = gRNAs, genomeSeqFile = genomeSeqFile), searchHits_args))
  }
    ## 3-2: With bulge:
  hits_bulge <- data.frame()
  if (findOffTargetsWithBulge) {
    ### Step3-1: search for off targets with bulge: start ###
    inform("Searching for off-targets with bulges ...")
    df_bulge <- do.call(getOfftargetWithBulge, c(list(gRNA_PAM = gRNAs), getOfftargetWithBulge_args))
    inform("Searching for off-targets with bulges: done!")
    
    ### Step3-2: filter off targets with bulge: start ###
    inform("Filtering off-targets with bulges ...")
    df_bulge <- do.call(filterCasOffinder, c(list(df_bulge = df_bulge), filterCasOffinder_args))
    inform("Filtering off-targets with bulges: done!")
    
    ### Step3-3: convert into hits-like object
    # strand, chrom, chromStart, chromEnd, name, gRNAplusPAM, gRNAPlusPAM_bulge, OffTargetSequence_bulge, n.mismatch, n.DNABulge, n.RNABulge, forViewInUCSC, score
    inform("Formatting off-targets with bulges: ...")
    hits_bulge <- formatCasOffinder(df_bulge, gRNAs = gRNAs)
    inform("Formatting off-targets with bulges: done!")
  }
  inform("Searching for off-targets: done!")
  ### Step3: search for off targets: finish ###
  
  ### Step4: filter offtargets: start ###
  if (dim(hits)[1] == 0 && dim(hits_bulge)[1] == 0) {
    summary <- do.call(getgRNASummary, c(list(gRNAs = gRNAs), getgRNASummary_args))
    return(summary)
  } else {
    ## 4-1a: get offtarget scores
    scores <- do.call(getOfftargetScoreWrap, c(list(hits = hits), getOfftargetScoreWrap_args))
    ## 4-1b: get offtarget scores for offtargets with bulge
    scores_bulge <- do.call(getOfftargetScoreBulgeWrap, c(list(hits = hits_bulge), getOfftargetScoreWrap_args))
    
    ## 4-2: annotate offtargets 
    scores <- combineScoresBulge(scores = scores, scores_bulge = scores_bulge)
    tem_res <- do.call(getOfftargetSummary, c(list(scores = scores), getOfftargetSummary_args))
    offTargets <- tem_res[[1]] # OfftargetAnalysis.xlsx
    summary <- tem_res[[2]] # Summary.xlsx
    on.target <- getOntarget(offTargets)

    ## 4-3: add paired and RE info to summary
    summary <- addInfoSummary(findgRNAs = findgRNAs,
                              findPairedgRNAOnly = findPairedgRNAOnly,
                              annotatePaired = annotatePaired,
                              summary = summary,
                              pairedInformation = pairedInformation,
                              REcutDetails = REcutDetails)
    
    ## 4-4: write gRNAs to bed file 
    gRNA.bed <- writegRNAsToBed(summary = summary,
                                offTargets = offTargets,
                                bedFile = bedFile,
                                useScore = useScore,
                                overlap.gRNA.positions = overlap.gRNA.positions)
    
    ## 4-5: update summary: add on.target, REs.isUnique100/50
    summary <- updateSummary(summary = summary, 
                             outputDir = outputDir,
                             offTargets = offTargets,
                             useEfficacyFromInputSeq = useEfficacyFromInputSeq,
                             efficacyFile = efficacyFile)
    
    ## 4-6: update summary: scan for REsites in flanking seq
    summary <- scanREsites(summary = summary, 
                           REcutDetails = REcutDetails, 
                           offTargets = offTargets, 
                           BSgenomeName = BSgenomeName, 
                           outputUniqueREs = outputUniqueREs)
    
    ## 4-7: update summary by adding free energy calculation
    summary <- updateSummaryFoldgRNAs(foldgRNAs = foldgRNAs, 
                                      summary = summary, 
                                      gRNA.size = gRNA.size, 
                                      gRNA.backbone = gRNA.backbone, 
                                      temperature = temperature)
    
    ## 4-8: update summary by adding gRNAs with no target in genome
    summary <- updateSummaryNoTargetgRNAs(summary = summary,
                                          gRNAs = gRNAs, 
                                          gRNA.size = gRNA.size,
                                          PAM = PAM,
                                          PAM.size = PAM.size,
                                          PAM.location = PAM.location,
                                          on.target = on.target,
                                          outputDir = outputDir)
    ## 4-9: add Lindel results 
    tem_res <- addLindelRes(predIndelFreq = predIndelFreq,
                            predictIndelFreq.onTargetOnly = predictIndelFreq.onTargetOnly,
                            offTargets = offTargets,
                            on.target = on.target,
                            BSgenomeName = BSgenomeName,
                            baseBeforegRNA.indelFreq = baseBeforegRNA.indelFreq,
                            baseAfterPAM.indelFreq = baseAfterPAM.indelFreq,
                            method.indelFreq = method.indelFreq,
                            summary = summary,
                            outputDir = outputDir)
    summary <- tem_res[[1]]
    indelFreq <- tem_res[[2]]
    fs2 <- tem_res[[3]]
  } 
  ### Step4: filter offtargets: finish ###
  
  ### Step5: return a list of results: start ###
  inform(paste0("All done! Please check output files in the directory: ", outputDir))
  list(on.target = on.target, 
       summary = summary, 
       offtarget = offTargets$offtargets,
       gRNAs.bedFormat = gRNA.bed, 
       REcutDetails = REcutDetails,
       REs.isUnique100 = summary$REs.isUnique100, 
       REs.isUnique50 = summary$REs.isUnique50,
       indelFreq = indelFreq, 
       fs2 = fs2)
  ### Step5: return a list of results: finish ###
}

