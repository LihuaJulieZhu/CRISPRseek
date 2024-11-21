test_that("test_offTargetAnalysis", {
  outputDir <- tempdir()
  chroms <- c("chr22", "chrX", "chr15", "chr18", "chr20", "chr5", "chr17", "chr19", "chr16", "chr10", "chr9", "chr1", "chr8", "chr12", "chrY", "chr2", "chr14", "chr11", "chr3", "chr7", "chr13", "chr4", "chr6")

  test.gRNAPlusPAM <- FALSE
  
  inputFilePath <- system.file("extdata", "inputseq.fa", package = "CRISPRseek")
  REpatternFile <- system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek")
  gRNAFilePath <- system.file("extdata", "testHsap_GATA1_ex2_gRNA1.fa", package = "CRISPRseek")
  
  # TEST1: paired search with RE:
  REcutDetails.pairedSearchRE <- read.xlsx(system.file("extdata/pairedSearch/withRE/REcutDetails.xlsx", package = "CRISPRseek"))
  pairedgRNAs.pairedSearchRE <- read.xlsx(system.file("extdata/pairedSearch/withRE/pairedgRNAs.xlsx", package = "CRISPRseek"))
  summary.pairedSearchRE <- read.xlsx(system.file("extdata/pairedSearch/withRE/Summary.xlsx", package = "CRISPRseek"))
  offtarget.pairedSearchRE  <- read.xlsx(system.file("extdata/pairedSearch/withRE/OfftargetAnalysis.xlsx", package = "CRISPRseek"))
  
  # TEST2: paired search without RE:
  REcutDetails.pairedSearchNoRE <- read.xlsx(system.file("extdata/pairedSearch/withoutRE/REcutDetails.xlsx", package = "CRISPRseek"))
  pairedgRNAs.pairedSearchNoRE <- read.xlsx(system.file("extdata/pairedSearch/withoutRE/pairedgRNAs.xlsx", package = "CRISPRseek"))
  summary.pairedSearchNoRE <- read.xlsx(system.file("extdata/pairedSearch/withoutRE/Summary.xlsx", package = "CRISPRseek"))
  offtarget.pairedSearchNoRE  <- read.xlsx(system.file("extdata/pairedSearch/withoutRE/OfftargetAnalysis.xlsx", package = "CRISPRseek"))
  
  # Unpaired search with RE:
  REcutDetails.unPairedSearchRE <- read.xlsx(system.file("extdata/unPairedSearch/withRE/REcutDetails.xlsx", package = "CRISPRseek"))
  pairedgRNAs.unPairedSearchRE <- read.xlsx(system.file("extdata/unPairedSearch/withRE/pairedgRNAs.xlsx", package = "CRISPRseek"))
  summary.unPairedSearchRE <- read.xlsx(system.file("extdata/unPairedSearch/withRE/Summary.xlsx", package = "CRISPRseek"))
  offtarget.unPairedSearchRE  <- read.xlsx(system.file("extdata/unPairedSearch/withRE/OfftargetAnalysis.xlsx", package = "CRISPRseek"))
  
  # Unpaired search without RE:
  REcutDetails.unPairedSearchNoRE <- read.xlsx(system.file("extdata/unPairedSearch/withoutRE/REcutDetails.xlsx", package = "CRISPRseek"))
  pairedgRNAs.unPairedSearchNoRE <- read.xlsx(system.file("extdata/unPairedSearch/withoutRE/pairedgRNAs.xlsx", package = "CRISPRseek"))
  summary.unPairedSearchNoRE <- read.xlsx(system.file("extdata/unPairedSearch/withoutRE/Summary.xlsx", package = "CRISPRseek"))
  offtarget.unPairedSearchNoRE  <- read.xlsx(system.file("extdata/unPairedSearch/withoutRE/OfftargetAnalysis.xlsx", package = "CRISPRseek"))
  
  # gRNAs provided search:
  REcutDetails.gRNAProvided <- read.xlsx(system.file("extdata/gRNAProvidedSearch/REcutDetails.xlsx", package = "CRISPRseek"))
  summary.gRNAProvided <- read.xlsx(system.file("extdata/gRNAProvidedSearch/Summary.xlsx", package = "CRISPRseek"))
  offtarget.gRNAProvided  <- read.xlsx(system.file("extdata/gRNAProvidedSearch/OfftargetAnalysis.xlsx", package = "CRISPRseek"))

  summary.pairedSearchRE <- summary.pairedSearchRE[, 1:22]
  summary.pairedSearchNoRE <- summary.pairedSearchNoRE[, 1:22]
  summary.unPairedSearchRE <- summary.unPairedSearchRE[, 1:21]
  summary.unPairedSearchNoRE <- summary.unPairedSearchNoRE[, 1:21]
  summary.gRNAProvided <- summary.gRNAProvided[, 1:21]

  if (!test.gRNAPlusPAM) {
    exclude.sum.col <-  grep("gRNAsPlusPAM", colnames(summary.gRNAProvided))
    exclude.oft.col <- grep("gRNAPlusPAM", colnames(offtarget.gRNAProvided))
    
    summary.pairedSearchRE <- summary.pairedSearchRE[, -exclude.sum.col]
    offtarget.pairedSearchRE <- offtarget.pairedSearchRE[, -exclude.oft.col]
    summary.pairedSearchNoRE <- summary.pairedSearchNoRE[, -exclude.sum.col]
    offtarget.pairedSearchNoRE <- offtarget.pairedSearchNoRE[, -exclude.oft.col]
    
    offtarget.unPairedSearchRE <- offtarget.unPairedSearchRE[, -exclude.oft.col]
    summary.unPairedSearchRE <- summary.unPairedSearchRE[, -exclude.sum.col]
    offtarget.unPairedSearchNoRE <- offtarget.unPairedSearchNoRE[, -exclude.oft.col]
    summary.unPairedSearchNoRE <- summary.unPairedSearchNoRE[, -exclude.sum.col]
    
    offtarget.gRNAProvided <- offtarget.gRNAProvided[, -exclude.oft.col]
    summary.gRNAProvided <- summary.gRNAProvided[, -exclude.sum.col]
  }
    
  # Actual testing:
  REcutDetails <- NULL
  pairedgRNAs <- NULL
  summary <- NULL
  offtarget <- NULL
  
  helper <- function(isRE = NULL, isPaired = NULL) {
    offTargetAnalysis(inputFilePath = inputFilePath, 
                      findgRNAsWithREcutOnly = isRE, 
                      findPairedgRNAOnly = isPaired,
                      REpatternFile = REpatternFile,
                      BSgenomeName = Hsapiens, 
                      txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                      chromToSearch = chroms,
                      orgAnn = org.Hs.egSYMBOL, 
                      max.mismatch = 3, 
                      min.score = 0.5, 
                      topN = 100,
                      outputDir= outputDir, 
                      overwrite = TRUE)
    REcutDetails <<- read.xlsx(file.path(outputDir, "REcutDetails.xlsx"))
    pairedgRNAs <<- read.xlsx(file.path(outputDir, "pairedgRNAs.xlsx"))
    summary <<- read.xlsx(file.path(outputDir, "Summary.xlsx"))
    offtarget <<- read.xlsx(file.path(outputDir, "OfftargetAnalysis.xlsx"))

    if (!test.gRNAPlusPAM) {
      summary <<- summary[, -exclude.sum.col]
      offtarget <<- offtarget[, -exclude.oft.col]
    }
  }
  
  # TEST1:
  inform("Testing for paired with RE sites ...")
  helper(isPaired = TRUE, isRE = TRUE)
  summary.pairedSearchRE$REname <- sapply(summary.pairedSearchRE$REname, function(x) ifelse(is.na(x), "", x))
  summary.pairedSearchRE$uniqREin100 <- sapply(summary.pairedSearchRE$uniqREin100, function(x) ifelse(is.na(x), "", x))
  summary.pairedSearchRE$uniqREin200 <- sapply(summary.pairedSearchRE$uniqREin200, function(x) ifelse(is.na(x), "", x))
  offtarget.pairedSearchRE$symbol <- sapply(offtarget.pairedSearchRE$symbol, function(x) ifelse(is.na(x), "", x))
  offtarget.pairedSearchRE$mismatch.distance2PAM <- sapply(offtarget.pairedSearchRE$mismatch.distance2PAM, function(x) ifelse(is.na(x), "", x)) # read.xlsx somehow convert "" into NA.
  
  expect_equal(REcutDetails.pairedSearchRE, REcutDetails)
  expect_equal(pairedgRNAs.pairedSearchRE, pairedgRNAs)
  expect_equal(summary.pairedSearchRE, summary[names(summary.pairedSearchRE)], tolerance = 0.01)
  expect_equal(offtarget.pairedSearchRE, offtarget[names(offtarget.pairedSearchRE)], tolerance = 0.001)
  
  # TEST2:
  inform("Testing for paired without RE sites ...")
  helper(isPaired = TRUE, isRE = FALSE)
  summary.pairedSearchNoRE$REname <- sapply(summary.pairedSearchNoRE$REname, function(x) ifelse(is.na(x), "", x))
  summary.pairedSearchNoRE$uniqREin100 <- sapply(summary.pairedSearchNoRE$uniqREin100, function(x) ifelse(is.na(x), "", x))
  summary.pairedSearchNoRE$uniqREin200 <- sapply(summary.pairedSearchNoRE$uniqREin200, function(x) ifelse(is.na(x), "", x))
  offtarget.pairedSearchNoRE$symbol <- sapply(offtarget.pairedSearchNoRE$symbol, function(x) ifelse(is.na(x), "", x))
  offtarget.pairedSearchNoRE$mismatch.distance2PAM <- sapply(offtarget.pairedSearchNoRE$mismatch.distance2PAM, function(x) ifelse(is.na(x), "", x)) # read.xlsx somehow convert "" into NA.
  
  expect_equal(REcutDetails.pairedSearchNoRE, REcutDetails)
  expect_equal(pairedgRNAs.pairedSearchNoRE, pairedgRNAs)
  expect_equal(summary.pairedSearchNoRE, summary[names(summary.pairedSearchNoRE)], tolerance = 0.01)
  expect_equal(offtarget.pairedSearchNoRE, offtarget[names(offtarget.pairedSearchNoRE)], tolerance = 0.001)
  
  # TEST3:
  inform("Testing for unpaired with RE sites ...")
  helper(isPaired = FALSE, isRE = TRUE)
  summary.unPairedSearchRE$REname <- sapply(summary.unPairedSearchRE$REname, function(x) ifelse(is.na(x), "", x))
  summary.unPairedSearchRE$uniqREin100 <- sapply(summary.unPairedSearchRE$uniqREin100, function(x) ifelse(is.na(x), "", x))
  summary.unPairedSearchRE$uniqREin200 <- sapply(summary.unPairedSearchRE$uniqREin200, function(x) ifelse(is.na(x), "", x))
  offtarget.unPairedSearchRE$symbol <- sapply(offtarget.unPairedSearchRE$symbol, function(x) ifelse(is.na(x), "", x))
  offtarget.unPairedSearchRE$mismatch.distance2PAM <- sapply(offtarget.unPairedSearchRE$mismatch.distance2PAM, function(x) ifelse(is.na(x), "", x)) # read.xlsx somehow convert "" into NA.
  
  expect_equal(REcutDetails.unPairedSearchRE, REcutDetails)
  expect_equal(pairedgRNAs.unPairedSearchRE, pairedgRNAs)
  expect_equal(summary.unPairedSearchRE, summary[names(summary.unPairedSearchRE)], tolerance = 0.01)
  expect_equal(offtarget.unPairedSearchRE, offtarget[names(offtarget.unPairedSearchRE)], tolerance = 0.001)
  
  # TEST4:
  inform("Testing for unpaired without RE sites ...")
  helper(isPaired = FALSE, isRE = FALSE)
  summary.unPairedSearchNoRE$REname <- sapply(summary.unPairedSearchNoRE$REname, function(x) ifelse(is.na(x), "", x))
  summary.unPairedSearchNoRE$uniqREin100 <- sapply(summary.unPairedSearchNoRE$uniqREin100, function(x) ifelse(is.na(x), "", x))
  summary.unPairedSearchNoRE$uniqREin200 <- sapply(summary.unPairedSearchNoRE$uniqREin200, function(x) ifelse(is.na(x), "", x))
  offtarget.unPairedSearchNoRE$symbol <- sapply(offtarget.unPairedSearchNoRE$symbol, function(x) ifelse(is.na(x), "", x))
  offtarget.unPairedSearchNoRE$mismatch.distance2PAM <- sapply(offtarget.unPairedSearchNoRE$mismatch.distance2PAM, function(x) ifelse(is.na(x), "", x)) # read.xlsx somehow convert "" into NA.
  
  expect_equal(REcutDetails.unPairedSearchNoRE, REcutDetails)
  expect_equal(pairedgRNAs.unPairedSearchNoRE, pairedgRNAs)
  expect_equal(summary.unPairedSearchNoRE, summary[names(summary.unPairedSearchNoRE)], tolerance = 0.01)
  expect_equal(offtarget.unPairedSearchNoRE, offtarget[names(offtarget.unPairedSearchNoRE)], tolerance = 0.001)
  
  # TEST5:
  inform("Testing for findgRNAs = FALSE ...")
  offTargetAnalysis(inputFilePath = gRNAFilePath, 
                    findgRNAs = FALSE, 
                    findgRNAsWithREcutOnly = FALSE, 
                    REpatternFile = REpatternFile, 
                    findPairedgRNAOnly = FALSE, 
                    BSgenomeName = Hsapiens, 
                    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                    orgAnn = org.Hs.egSYMBOL, 
                    max.mismatch = 3, 
                    chromToSearch = c("chrX", "chr11", "chr4", "chr5", "chr6", "chr8", "chr15", "chr3"),
                    min.score = 0.5, 
                    topN = 100,
                    outputDir = outputDir, 
                    overwrite = TRUE)
    
  REcutDetails <- read.xlsx(file.path(outputDir, "REcutDetails.xlsx"))
  summary <- read.xlsx(file.path(outputDir, "Summary.xlsx"))
  offtarget <- read.xlsx(file.path(outputDir, "OfftargetAnalysis.xlsx"))
  if (!test.gRNAPlusPAM) {
    summary <- summary[, -exclude.sum.col]
    offtarget <- offtarget[, -exclude.oft.col]
  }

  expect_equal(REcutDetails.gRNAProvided, REcutDetails)
  
  summary.gRNAProvided$REname <- sapply(summary.gRNAProvided$REname, function(x) ifelse(is.na(x), "", x))
  summary.gRNAProvided$uniqREin100 <- sapply(summary.gRNAProvided$uniqREin100, function(x) ifelse(is.na(x), "", x))
  summary.gRNAProvided$uniqREin200 <- sapply(summary.gRNAProvided$uniqREin200, function(x) ifelse(is.na(x), "", x))
  offtarget.gRNAProvided$symbol <- sapply(offtarget.gRNAProvided$symbol, function(x) ifelse(is.na(x), "", x))
  offtarget.gRNAProvided$mismatch.distance2PAM <- sapply(offtarget.gRNAProvided$mismatch.distance2PAM, function(x) ifelse(is.na(x), "", x)) # read.xlsx somehow convert "" into NA.
  
  expect_equal(summary.gRNAProvided, summary[names(summary.gRNAProvided)], tolerance = 0.01)
  expect_equal(offtarget.gRNAProvided, offtarget[names(offtarget.gRNAProvided)], tolerance = 0.001)

})
