test_that("test_targetOffNotFound5primePAM", {
  cat(paste0(
    "Testing for 5 prime PAM target and offtarget not found when q",
    "gRNA sequence is from another species serving as a negative control…\n"
  ))
  
  inputSeq <- DNAStringSet("CAGTATTACTGATATTGGTTTTTCCAGGG")
  outputDir <- getwd()
  BSgenomeName <- Mmusculus
  chroms <- "chr12"

  expect_warning(offTargetAnalysis(inputSeq,
    gRNAoutputName = "negativeCntl",
    PAM = "GTN",
    PAM.pattern = "^GTN",
    PAM.location = "5prime",
    format = "fasta",
    findgRNAs = TRUE,
    findgRNAsWithREcutOnly = FALSE,
    findPairedgRNAOnly = FALSE,
    annotatePaired = FALSE,
    BSgenomeName = BSgenomeName,
    chromToSearch = chroms, max.mismatch = 0, annotateExon = FALSE,
    fetchSequence = FALSE, outputDir = outputDir, overwrite = TRUE
  ))
  summary <- read.table("Summary.xls",
    sep = "\t", header = TRUE,
    stringsAsFactors = FALSE
  )

  expect_equal(
    "perfect match not found",
    summary$top1Hit.onTarget.MMdistance2PAM
  )
  expect_equal(NA, summary$top5OfftargetTotalScore)
  expect_equal("GTNTTACTGATATTGGTTTTTCC", summary$gRNAsPlusPAM)
  
  cat("2nd test")
  inputSeq <- DNAStringSet("GTATTACTGATATTGGTGGGAGG")

  offTargetAnalysis(inputSeq,
                                   gRNAoutputName = "negativeCntl",
                                   PAM = "GTN",
                                   PAM.pattern = "^NTN",
                                   PAM.location = "5prime",
                                   format = "fasta",
                                   findgRNAs = TRUE,
                                   findgRNAsWithREcutOnly = FALSE,
                                   findPairedgRNAOnly = FALSE,
                                   annotatePaired = FALSE,
                                   BSgenomeName = BSgenomeName,
                                   chromToSearch = chroms, 
                                   max.mismatch = 4, 
                                   annotateExon = FALSE,
                                   fetchSequence = FALSE, 
                                   outputDir = outputDir, 
                                   overwrite = TRUE
  )
  summary <- read.table("Summary.xls",
                        sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE
  )
  
  expect_equal(
    "perfect match not found",
    summary$top1Hit.onTarget.MMdistance2PAM
  )
  expect_equal(3.6, summary$top5OfftargetTotalScore, tolerance = 0.1)
  expect_equal("GTNTTACTGATATTGGTGGGAGG", summary$gRNAsPlusPAM)
  
  
})
