test_that("test_mouse_orgAnn_annotateExon",  {
  inputFilePath = system.file("extdata", "RIPK1stop.fa", package = "CRISPRseek")
  REpatternFile <- system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek")
  
  test.gRNAPlusPAM <- TRUE
  outputDir <- tempdir()

  summary.NOorgAnn <- read.xlsx(system.file("extdata/testMouse/NOorgAnn/Summary.xlsx", package = "CRISPRseek"))
  offtarget.NOorgAnn  <- read.xlsx(system.file("extdata/testMouse/NOorgAnn/OfftargetAnalysis.xlsx", package = "CRISPRseek"))
  summary.WithorgAnn <- read.xlsx(system.file("extdata/testMouse/WithorgAnn/Summary.xlsx", package = "CRISPRseek"))
  offtarget.WithorgAnn  <- read.xlsx(system.file("extdata/testMouse/WithorgAnn/OfftargetAnalysis.xlsx", package = "CRISPRseek"))
  summary.NOTannotateExon <- read.xlsx(system.file("extdata/testMouse/NOTannotateExon/Summary.xlsx", package = "CRISPRseek"))
  offtarget.NOTannotateExon  <- read.xlsx(system.file("extdata/testMouse/NOTannotateExon/OfftargetAnalysis.xlsx", package = "CRISPRseek"))

  summary.WithorgAnn <- summary.WithorgAnn[, 1:21]
  summary.NOorgAnn <- summary.NOorgAnn [, 1:21]
  summary.NOTannotateExon <- summary.NOTannotateExon[, 1:21]

  if (!test.gRNAPlusPAM) {
    exclude.sum.col <-  grep("gRNAsPlusPAM", colnames(summary.WithorgAnn))
    exclude.oft.col <- grep("gRNAPlusPAM", colnames(offtarget.WithorgAnn)) 
    offtarget.NOTannotateExon <- offtarget.NOTannotateExon[, -exclude.oft.col]
    summary.NOTannotateExon <- summary.NOTannotateExon[, -exclude.sum.col]
    offtarget.WithorgAnn <- offtarget.WithorgAnn[, -exclude.oft.col]
    summary.WithorgAnn <- summary.WithorgAnn[, -exclude.sum.col]
    summary.NOorgAnn <- summary.NOorgAnn[, -exclude.sum.col]
    offtarget.NOorgAnn <- offtarget.NOorgAnn[, -exclude.oft.col]
  }
  chroms <- c("chr1", "chr5", "chr13")
  
  # TEST1: mouse with orgAnn
  inform("Testing for mouse with orgAnn ...")
  offTargetAnalysis(inputFilePath, 
                    findgRNAs = TRUE, 
                    findgRNAsWithREcutOnly = FALSE, 
                    findPairedgRNAOnly = FALSE,
                    BSgenomeName = Mmusculus, 
                    annotateExon = TRUE, 
                    chromToSearch = chroms,
                    min.score = 0.5, 
                    topN = 100, 
                    outputDir = outputDir, 
                    overwrite = TRUE, 
                    max.mismatch = 1, 
                    txdb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                    orgAnn = org.Mm.egSYMBOL)
  
  summary <- read.xlsx(file.path(outputDir, "Summary.xlsx"))
  offtarget <- read.xlsx(file.path(outputDir, "OfftargetAnalysis.xlsx"))

  if (!test.gRNAPlusPAM) {
    summary <- summary[, -exclude.sum.col]
    offtarget <- offtarget[, -exclude.oft.col]
  }
    
  summary.WithorgAnn$REname <- sapply(summary.WithorgAnn$REname, function(x) ifelse(is.na(x), "", x))
  summary.WithorgAnn$uniqREin200 <- sapply(summary.WithorgAnn$uniqREin200, function(x) ifelse(is.na(x), "", x))
  summary.WithorgAnn$uniqREin100 <- sapply(summary.WithorgAnn$uniqREin100, function(x) ifelse(is.na(x), "", x))
  offtarget.WithorgAnn$symbol <- sapply(offtarget.WithorgAnn$symbol, function(x) ifelse(is.na(x), "", x))
  for (x in paste("topOfftarget", 1:10, "MMdistance2PAM", sep = "")) {
    summary.WithorgAnn[, x] <- as.character(summary.WithorgAnn[, x])
    summary[, x] <- as.character(summary[, x])
  }
  offtarget.WithorgAnn$mismatch.distance2PAM <- sapply(offtarget.WithorgAnn$mismatch.distance2PAM, function(x) ifelse(is.na(x), "", x))
    
  expect_equal(summary.WithorgAnn, summary[names(summary.WithorgAnn)], tolerance = 0.01)
  expect_equal(offtarget.WithorgAnn, offtarget[names(offtarget.WithorgAnn)], tolerance = 0.01)
  
  # TEST2: mouse without orgAnn
  inform("Testing for mouse without orgAnn ...")
  expect_warning(offTargetAnalysis(inputFilePath, 
                                   findgRNAs = TRUE, 
                                   findgRNAsWithREcutOnly = FALSE, 
                                   findPairedgRNAOnly = FALSE,
                                   BSgenomeName = Mmusculus, 
                                   annotateExon = TRUE, 
                                   chromToSearch = chroms,
                                   min.score = 0.5, 
                                   topN = 100, 
                                   outputDir = outputDir, 
                                   overwrite = TRUE, 
                                   max.mismatch = 1, 
                                   txdb = TxDb.Mmusculus.UCSC.mm10.knownGene))

  summary <- read.xlsx(file.path(outputDir, "Summary.xlsx"))
  offtarget <- read.xlsx(file.path(outputDir, "OfftargetAnalysis.xlsx"))
  if (!test.gRNAPlusPAM) {
    summary <- summary[, -exclude.sum.col]
    offtarget <- offtarget[, -exclude.oft.col]
  }
  
  summary.NOorgAnn$REname <- sapply(summary.NOorgAnn$REname, function(x) ifelse(is.na(x), "", x))
  summary.NOorgAnn$uniqREin200 <- sapply(summary.NOorgAnn$uniqREin200, function(x) ifelse(is.na(x), "", x))
  summary.NOorgAnn$uniqREin100 <- sapply(summary.NOorgAnn$uniqREin100, function(x) ifelse(is.na(x), "", x))
  for (x in paste("topOfftarget", 1:10, "MMdistance2PAM", sep = "")) {
    summary.NOorgAnn[, x] <- as.character(summary.NOorgAnn[, x])
    summary[, x] <- as.character(summary[, x])
  }
  offtarget.NOorgAnn$mismatch.distance2PAM <- sapply(offtarget.NOorgAnn$mismatch.distance2PAM, function(x) ifelse(is.na(x), "", x))
  
  expect_equal(summary.NOorgAnn, summary[names(summary.NOorgAnn)], tolerance = 0.01)
  expect_equal(offtarget.NOorgAnn, offtarget[names(offtarget.NOorgAnn)], tolerance = 0.01)
  
  # TEST3: mouse without orgAnn and without annotateExon
  inform("Testing for mouse without orgAnn ...")
  offTargetAnalysis(inputFilePath, 
                    findgRNAs = TRUE, 
                    findgRNAsWithREcutOnly = FALSE, 
                    findPairedgRNAOnly = FALSE,
                    BSgenomeName = Mmusculus, 
                    annotateExon = FALSE, 
                    chromToSearch = chroms,
                    min.score = 0.5, 
                    topN = 100, 
                    outputDir = outputDir, 
                    overwrite = TRUE, 
                    max.mismatch = 1)

  summary <- read.xlsx(file.path(outputDir, "Summary.xlsx"))
  offtarget <- read.xlsx(file.path(outputDir, "OfftargetAnalysis.xlsx"))
  if (!test.gRNAPlusPAM) {
    summary <- summary[, -exclude.sum.col]
    offtarget <- offtarget[, -exclude.oft.col]
  }
  
  summary.NOTannotateExon$REname <- sapply(summary.NOTannotateExon$REname, function(x) ifelse(is.na(x), "", x))
  summary.NOTannotateExon$uniqREin100 <- sapply(summary.NOTannotateExon$uniqREin100, function(x) ifelse(is.na(x), "", x))
  summary.NOTannotateExon$uniqREin200 <- sapply(summary.NOTannotateExon$uniqREin200, function(x) ifelse(is.na(x), "", x))
  for (x in paste("topOfftarget", 1:10, "MMdistance2PAM", sep = "")) {
    summary.NOTannotateExon[, x] <- as.character(summary.NOTannotateExon[, x])
    summary[, x] <- as.character(summary[, x])
  }
  offtarget.NOTannotateExon$mismatch.distance2PAM <- sapply(offtarget.NOTannotateExon$mismatch.distance2PAM, function(x) ifelse(is.na(x), "", x))
  
  expect_equal(summary.NOTannotateExon, summary[names(summary.NOTannotateExon)], tolerance = 0.01)
  expect_equal(offtarget.NOTannotateExon, offtarget[names(offtarget.NOTannotateExon)], tolerance = 0.01)
})
