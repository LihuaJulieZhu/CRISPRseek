inputFilePath <- system.file("extdata/testNoTarget/",
    "negativeCntl-gRNA.fa", package = "CRISPRseek")

outputDir <- getwd()

library(CRISPRseek)
library("BSgenome.Mmusculus.UCSC.mm10")

BSgenomeName = Mmusculus


summary.truth <- read.table(
    system.file("extdata/testNoTarget/",
    "Summary.xls", package = "CRISPRseek"), sep = "\t", header = TRUE,
    stringsAsFactors = FALSE)
offtarget.truth  <- read.table(
    system.file("extdata/testNoTarget",
    "OfftargetAnalysis.xls", package = "CRISPRseek"), sep = "\t",
    header = TRUE, stringsAsFactors = FALSE)



chroms <- c("chr19", "chr18", "chr13",  "chr12", "chr11")

test_targetNotFound <- function() {
     cat("Testing for target not found when a gRNA sequence is from another species serving as a negative control…\n")

     offTargetAnalysis(inputFilePath,
                   scoring.method = "CFDscore",
                   min.score = 0.001,
                   format = "fasta",
                   findgRNAs = FALSE, 
                   findgRNAsWithREcutOnly = FALSE, 
                   findPairedgRNAOnly = FALSE, 
                   annotatePaired = FALSE,
                   BSgenomeName = BSgenomeName,
                   chromToSearch= chroms, 
                   max.mismatch = 3,
                   annotateExon = FALSE,
                   fetchSequence = FALSE, 
                   outputDir = outputDir,
                   overwrite = TRUE)

    summary <- read.table("Summary.xls", sep = "\t", header = TRUE,
        stringsAsFactors = FALSE)
    offtarget <- read.table("OfftargetAnalysis.xls", sep = "\t", header = TRUE,
        stringsAsFactors = FALSE)

    if (checkEquals(summary.truth, summary, tolerance = 0.01))
        cat("Summary passed test for test_targetNotFound\n")
    else
        cat("Summary failed test for gRNA with no on-target found\n")
    if (checkEquals(offtarget.truth, offtarget,tolerance = 0.001))
        cat("off target analysis details passed test for test_targetNotFound\n")
    else
        cat("off target analysis details failed test for gRNA with no on-target found\n")
}
