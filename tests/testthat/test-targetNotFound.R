test_that("test_targetNotFound", {
    cat("Testing for target not found when a gRNA sequence is from another species serving as a negative control…\n")

  inputFilePath <- system.file("extdata/testNoTarget/",
                               "negativeCntl-gRNA.fa", package = "CRISPRseek")
  outputDir <- getwd()
  
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

     offTargetAnalysis(inputFilePath, scoring.method = "CFDscore", 
        min.score = 0.001, format = "fasta", findgRNAs = FALSE, 
        findgRNAsWithREcutOnly = FALSE, findPairedgRNAOnly = FALSE, 
        annotatePaired = FALSE, BSgenomeName = BSgenomeName, 
        chromToSearch = chroms, max.mismatch = 3, annotateExon = FALSE, 
        fetchSequence = FALSE, outputDir = outputDir, overwrite = TRUE)
     summary <- read.table("Summary.xls", sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE)
     offtarget <- read.table("OfftargetAnalysis.xls", sep = "\t", 
        header = TRUE, stringsAsFactors = FALSE)
     expect_equal(summary.truth, summary, tolerance = 0.01)
     expect_equal(offtarget.truth, offtarget, tolerance = 0.001)
})
