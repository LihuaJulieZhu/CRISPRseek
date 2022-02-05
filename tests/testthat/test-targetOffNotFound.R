test_that("test_targetOffNotFound", {
     cat(paste0("Testing for target and offtarget not found when a gRNA",
       "sequence is from another species serving as a negative control…\n"))
     inputFilePath <- system.file("extdata/testNoTarget/",
                               "negativeCntl-gRNA.fa", package = "CRISPRseek")
     outputDir <- getwd()
     BSgenomeName = Mmusculus
     chroms <-  "chr12"

     expect_warning(offTargetAnalysis(inputFilePath, 
                                      scoring.method = "CFDscore",
        min.score = 0.001, format = "fasta", findgRNAs = FALSE,
        findgRNAsWithREcutOnly = FALSE, findPairedgRNAOnly = FALSE,
        annotatePaired = FALSE, BSgenomeName = BSgenomeName,
        chromToSearch = chroms, max.mismatch = 2, annotateExon = FALSE,
        fetchSequence = FALSE, outputDir = outputDir, overwrite = TRUE))
     summary <- read.table("Summary.xls", sep = "\t", header = TRUE,
        stringsAsFactors = FALSE)
     
     expect_equal("perfect match not found",
                  summary$top1Hit.onTarget.MMdistance2PAM)
     expect_equal(NA, summary$top5OfftargetTotalScore)
})
