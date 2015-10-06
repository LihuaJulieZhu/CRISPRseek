\name{CRISPRseek-package}
\alias{CRISPRseek-package}
\alias{CRISPRseek}
\docType{package}
\title{
Design of target-specific guide RNAs (gRNAs) in CRISPR-Cas9, genome-editing systems
}
\description{
Design of target-specific gRNAs for the CRISPR-Cas9 system by automatically 
finding potential gRNAs (paired/not paired), with/without restriction enzyme 
cut site(s) in a given sequence, searching for off targets with user defined 
maximum number of mismatches, calculating score of each off target based on 
mismatch positions in the off target and a penalty weight matrix, filtering off 
targets with user-defined criteria, and annotating off targets with flank 
sequences, whether located in exon or not. Summary report is also generated 
with gRNAs ranked by total topN off target score, annotated with restriction 
enzyme cut sites, gRNA efficacy and  possible paired gRNAs. Detailed paired gRNAs 
information and restriction enzyme cut sites are stored in separate files in 
the output directory specified by the user. In total, four tab delimited files 
are generated in the output directory: OfftargetAnalysis.xls (off target 
details), Summary.xls (gRNA summary), REcutDetails.xls (restriction enzyme 
cut sites of each gRNA), and pairedgRNAs.xls (potential paired gRNAs).
}
\details{
\tabular{ll}{
Package: \tab CRISPRseek\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-10-04\cr
License: \tab GPL (>= 2)\cr
}
Function offTargetAnalysis integrates all steps of off target analysis into one 
function call
}
\author{
Lihua Julie Zhu and Michael Brodsky
Maintainer: julie.zhu@umassmed.edu
}
\references{
Mali P, Aach J, Stranges PB, Esvelt KM, Moosburner M, Kosuri S, Yang L, 
Church GM.CAS9 transcriptional activators for target specificity screening 
and paired nickases for cooperative genome engineering. Nat Biotechnol. 2013. 
31(9):833-8
Patrick D Hsu, David A Scott, Joshua A Weinstein, F Ann Ran, Silvana Konermann, 
Vineeta Agarwala, Yinqing Li, Eli J Fine, Xuebing Wu, Ophir Shalem, Thomas J 
Cradick, Luciano A Marraffini, Gang Bao & Feng Zhang. DNA targeting specificity 
of rNA-guided Cas9 nucleases.  Nat Biotechnol. 2013. 31:827-834
Lihua Julie Zhu, Benjamin R. Holmes, Neil Aronin and Michael Brodsky. CRISPRseek: a Bioconductor package to identify target-specific guide RNAs for CRISPR-Cas9 genome-editing systems. Plos One Sept 23rd 2014
}

\keyword{ package }
\seealso{
offTargetAnalysis
}
\examples{
    library(CRISPRseek)
    library("BSgenome.Hsapiens.UCSC.hg19")
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
    outputDir <- getwd()
    inputFilePath <- system.file("extdata", "inputseq.fa", package = "CRISPRseek")
    REpatternFile <- system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek")
######## Scenario 1. Target and off-target analysis for paired gRNAs with 
######## one of the pairs overlap RE sites
    results <- offTargetAnalysis(inputFilePath, findgRNAsWithREcutOnly=TRUE,
        REpatternFile =REpatternFile,findPairedgRNAOnly=TRUE, 
        BSgenomeName=Hsapiens, txdb=TxDb.Hsapiens.UCSC.hg19.knownGene, 
        orgAnn = org.Hs.egSYMBOL,max.mismatch = 1, chromToSearch = "chrX",
        outputDir = outputDir,overwrite = TRUE)

######## Scenario 2. Target and off-target analysis for paired gRNAs with or 
######## without RE sites
    results <- offTargetAnalysis(inputFilePath, findgRNAsWithREcutOnly = FALSE,
        REpatternFile = REpatternFile,findPairedgRNAOnly = TRUE,
        BSgenomeName = Hsapiens, txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
        orgAnn = org.Hs.egSYMBOL,max.mismatch = 1, chromToSearch = "chrX",
        outputDir = outputDir, overwrite = TRUE)

######## Scenario 3. Target and off-target analysis for gRNAs overlap RE sites

    results <- offTargetAnalysis(inputFilePath, findgRNAsWithREcutOnly = TRUE,
        REpatternFile = REpatternFile,findPairedgRNAOnly = FALSE, 
        BSgenomeName = Hsapiens, txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
        orgAnn = org.Hs.egSYMBOL, max.mismatch = 1, chromToSearch = "chrX", 
        outputDir = outputDir, overwrite = TRUE)

######## Scenario 4. Off-target analysis for all potential gRNAs, this will 
########be the slowest among the aforementioned scenarios.

    results <- offTargetAnalysis(inputFilePath, findgRNAsWithREcutOnly = FALSE, 
        REpatternFile = REpatternFile,findPairedgRNAOnly = FALSE,
        BSgenomeName = Hsapiens, txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
         orgAnn = org.Hs.egSYMBOL, max.mismatch = 1, chromToSearch = "chrX", 
        outputDir = outputDir,overwrite = TRUE)

######## Scenario 5. Target and off-target analysis for gRNAs input by user.
    gRNAFilePath <- system.file("extdata", "testHsap_GATA1_ex2_gRNA1.fa",
        package="CRISPRseek")
    results <- offTargetAnalysis(inputFilePath = gRNAFilePath, findgRNAs = FALSE, 
        findgRNAsWithREcutOnly = FALSE, REpatternFile = REpatternFile,
        findPairedgRNAOnly = FALSE, BSgenomeName = Hsapiens, 
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
         orgAnn = org.Hs.egSYMBOL, max.mismatch = 1, chromToSearch = "chrX", 
        outputDir = outputDir, overwrite = TRUE)

####### Scenario 6. Quick gRNA finding  without target and off-target analysis 
    results <- offTargetAnalysis(inputFilePath, findgRNAsWithREcutOnly = TRUE,
        REpatternFile = REpatternFile,findPairedgRNAOnly = TRUE,  
        chromToSearch = "", outputDir = outputDir, overwrite = TRUE)

####### Scenario 7. Quick gRNA finding  with gRNA efficacy analysis
    results <- offTargetAnalysis(inputFilePath, findgRNAsWithREcutOnly = TRUE,
        REpatternFile = REpatternFile,findPairedgRNAOnly = TRUE,
	BSgenomeName = Hsapiens, annotateExon = FALSE,
        max.mismatch = 0, outputDir = outputDir, overwrite = TRUE)
}
