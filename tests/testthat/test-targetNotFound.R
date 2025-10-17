test_that("test_targetNotFound", {
  inform("Testing for target not found when a gRNA sequence is from another species serving as a negative control ...")

  inputFilePath <- system.file("extdata/testNoTarget/", "negativeCntl-gRNA.fa", package = "CRISPRseek")
  outputDir <- tempdir()
  BSgenomeName = Mmusculus
  
  summary.truth <- read.xlsx(system.file("extdata/testNoTarget/Summary.xlsx", package = "CRISPRseek"), na.strings = "")
  summary.truth$REs.isUnique100 <- NA_real_
  summary.truth$REs.isUnique50 <- NA_real_

  chroms <- c("chr19", "chr18", "chr13",  "chr12", "chr11")
  expect_warning(offTargetAnalysis(inputFilePath, 
                                   scoring.method = "CFDscore", 
                                   min.score = 0.001, 
                                   findgRNAs = FALSE, 
                                   findgRNAsWithREcutOnly = FALSE, 
                                   findPairedgRNAOnly = FALSE, 
                                   annotatePaired = FALSE, 
                                   BSgenomeName = BSgenomeName, 
                                   chromToSearch = chroms, 
                                   max.mismatch = 3, 
                                   annotateExon = FALSE, 
                                   fetchSequence = FALSE, 
                                   outputDir = outputDir, 
                                   overwrite = TRUE))

summary <- read.xlsx(file.path(outputDir, "Summary.xlsx"), na.strings = "")
expect_equal(summary.truth, summary[names(summary.truth)], tolerance = 0.01)

offtarget.truth  <- read.xlsx(system.file("extdata/testNoTarget/OfftargetAnalysis.xlsx", package = "CRISPRseek"), na.strings = "")
offtarget <- read.xlsx(file.path(outputDir, "OfftargetAnalysis.xlsx"), na.strings = "")

# Since 1.44.1, more columns regarding bulge analysis are being added to the xlsx:
expect_equal(offtarget.truth, offtarget[, names(offtarget.truth)], tolerance = 0.001)

})
