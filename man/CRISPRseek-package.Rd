\name{CRISPRseek-package}
\alias{CRISPRseek-package}
\alias{CRISPRseek}
\docType{package}
\title{
Seeking and evaluating target and off target sequences for the CRISPR-Cas9 
system
}
\description{
Target and off target analysis for the CRISPR-Cas9 system by automatically 
finding potential spacers (paired/not paired), with/without restriction enzyme 
cut site(s) in a given sequence, searching for off targets with user defined 
maximum number of mismatches, calculating score of each off target based on 
mismatch positions in the off target and a penalty weight matrix, filtering off 
targets with user-defined criteria, and annotating off targets with flank 
sequences, whether located in exon or not. Summary report is also generated 
with spacers ranked by total topN off target score, annotated with restriction 
enzyme cut sites and  possible paired spacers. Detailed paired spacers 
information and restriction enzyme cut sites are stored in separate files in 
the output directory specified by the user. In total, four tab delimited files 
are generated in the output directory: OfftargetAnalysis.xls (off target 
details), Summary.xls (spacer summary), REcutDetails.xls (restriction enzyme 
cut sites of each spacer), and pairedSpacers.xls (potential paired spacers).
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
}
\keyword{ package }
\seealso{
offTargetAnalysis
}
\examples{
    library(CRISPRseek)
    library("BSgenome.Hsapiens.UCSC.hg19")
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    outputDir <- getwd()
    inputFilePath <- system.file("extdata", "inputseq.fa", package = "CRISPRseek")
    REpatternFile <- system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek")
######## Scenario 1. Target and off-target analysis for paired spacers with 
######## one of the pairs overlap RE sites

    offTargetAnalysis(inputFilePath, findSpacersWithREcutOnly=TRUE,
        REpatternFile =REpatternFile,findPairedSpacerOnly=TRUE, 
        BSgenomeName=Hsapiens, txdb=TxDb.Hsapiens.UCSC.hg19.knownGene, 
        max.mismatch = 1, chromToSearch = "chrX",
        outputDir = outputDir,overwrite = TRUE)

######## Scenario 2. Target and off-target analysis for paired spacers with or 
######## without RE sites
    offTargetAnalysis(inputFilePath, findSpacersWithREcutOnly = FALSE,
        REpatternFile = REpatternFile,findPairedSpacerOnly = TRUE,
        BSgenomeName = Hsapiens, txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
        max.mismatch = 1, chromToSearch = "chrX",
        outputDir = outputDir, overwrite = TRUE)

######## Scenario 3. Target and off-target analysis for spacers overlap RE sites

    offTargetAnalysis(inputFilePath, findSpacersWithREcutOnly = TRUE,
        REpatternFile = REpatternFile,findPairedSpacerOnly = FALSE, 
        BSgenomeName = Hsapiens, txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
        max.mismatch = 1, chromToSearch = "chrX", 
        outputDir = outputDir, overwrite = TRUE)

######## Scenario 4. Off-target analysis for all potential spacers, this will 
########be the slowest among the aforementioned scenarios.

    offTargetAnalysis(inputFilePath, findSpacersWithREcutOnly = FALSE, 
        REpatternFile = REpatternFile,findPairedSpacerOnly = FALSE,
        BSgenomeName = Hsapiens, txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
        max.mismatch = 1, chromToSearch = "chrX", 
        outputDir = outputDir,overwrite = TRUE)

######## Scenario 5. Target and off-target analysis for spacers input by user.
    spacerFilePath <- system.file("extdata", "testHsap_GATA1_ex2_gRNA1.fa",
        package="CRISPRseek")
    offTargetAnalysis(inputFilePath = spacerFilePath, findSpacers = FALSE, 
        findSpacersWithREcutOnly = FALSE, REpatternFile = REpatternFile,
        findPairedSpacerOnly = FALSE, BSgenomeName = Hsapiens, 
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
        max.mismatch = 1, chromToSearch = "chrX", 
        outputDir = outputDir, overwrite = TRUE)

####### Scenario 6. Quick spacer finding  without target and off-target analysis 
    offTargetAnalysis(inputFilePath, findSpacersWithREcutOnly = TRUE,
        REpatternFile = REpatternFile,findPairedSpacerOnly = TRUE,  
        chromToSearch = "", outputDir = outputDir, overwrite = TRUE)
}
