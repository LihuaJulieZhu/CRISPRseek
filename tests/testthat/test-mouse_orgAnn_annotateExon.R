test_that("test_mouse_orgAnn_annotateExon",  {
  inputFilePath <- system.file("extdata", "RIPK1stop.fa", package = "CRISPRseek")
  REpatternFile <- system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek")
  
  test.gRNAPlusPAM <- TRUE
  outputDir <- tempdir()

  summary.WithorgAnn <- read.xlsx(system.file("extdata/testMouse/WithorgAnn/Summary.xlsx", package = "CRISPRseek"), na.strings = "")
  offtarget.WithorgAnn  <- read.xlsx(system.file("extdata/testMouse/WithorgAnn/OfftargetAnalysis.xlsx", package = "CRISPRseek"), na.strings = "")
  summary.NOorgAnn <- read.xlsx(system.file("extdata/testMouse/NOorgAnn/Summary.xlsx", package = "CRISPRseek"), na.strings = "")
  offtarget.NOorgAnn  <- read.xlsx(system.file("extdata/testMouse/NOorgAnn/OfftargetAnalysis.xlsx", package = "CRISPRseek"), na.strings = "")
  summary.NOTannotateExon <- read.xlsx(system.file("extdata/testMouse/NOTannotateExon/Summary.xlsx", package = "CRISPRseek"), na.strings = "")
  offtarget.NOTannotateExon  <- read.xlsx(system.file("extdata/testMouse/NOTannotateExon/OfftargetAnalysis.xlsx", package = "CRISPRseek"), na.strings = "")

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
  
  summary <- read.xlsx(file.path(outputDir, "Summary.xlsx"), na.strings = "")
  offtarget <- read.xlsx(file.path(outputDir, "OfftargetAnalysis.xlsx"), na.strings = "")

  if (!test.gRNAPlusPAM) {
    summary <- summary[, -exclude.sum.col]
    offtarget <- offtarget[, -exclude.oft.col]
  }
    
 suppressWarnings({
    summary.WithorgAnn$top5OfftargetTotalScore <- as.numeric(summary.WithorgAnn$top5OfftargetTotalScore)
    summary.WithorgAnn$top10OfftargetTotalScore <- as.numeric(summary.WithorgAnn$top10OfftargetTotalScore)
    summary$topOfftarget1MMdistance2PAM <- as.numeric(summary$topOfftarget1MMdistance2PAM)
    offtarget.WithorgAnn$mismatch.distance2PAM <- offtarget.WithorgAnn$mismatch.distance2PAM %>% as.numeric() %>% as.character()
  })
  
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

  summary <- read.xlsx(file.path(outputDir, "Summary.xlsx"), na.strings = "")
  offtarget <- read.xlsx(file.path(outputDir, "OfftargetAnalysis.xlsx"), na.strings = "")
  if (!test.gRNAPlusPAM) {
    summary <- summary[, -exclude.sum.col]
    offtarget <- offtarget[, -exclude.oft.col]
  }
  
  for (x in paste("topOfftarget", 1:10, "MMdistance2PAM", sep = "")) {
    summary.NOorgAnn[, x] <- as.character(summary.NOorgAnn[, x])
    summary[, x] <- as.character(summary[, x])
  }
  offtarget.NOorgAnn$mismatch.distance2PAM <- sapply(offtarget.NOorgAnn$mismatch.distance2PAM, function(x) ifelse(is.na(x), "", x))
  
  suppressWarnings({
    summary.NOorgAnn$top5OfftargetTotalScore <- as.numeric(summary.NOorgAnn$top5OfftargetTotalScore)
    summary.NOorgAnn$top10OfftargetTotalScore <- as.numeric(summary.NOorgAnn$top10OfftargetTotalScore)
    offtarget.NOorgAnn$mismatch.distance2PAM <- offtarget.NOorgAnn$mismatch.distance2PAM %>% as.numeric() %>% as.character()
  })
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

  summary <- read.xlsx(file.path(outputDir, "Summary.xlsx"), na.strings = "")
  offtarget <- read.xlsx(file.path(outputDir, "OfftargetAnalysis.xlsx"), na.strings = "")
  if (!test.gRNAPlusPAM) {
    summary <- summary[, -exclude.sum.col]
    offtarget <- offtarget[, -exclude.oft.col]
  }
  
  for (x in paste("topOfftarget", 1:10, "MMdistance2PAM", sep = "")) {
    summary.NOTannotateExon[, x] <- as.character(summary.NOTannotateExon[, x])
    summary[, x] <- as.character(summary[, x])
  }

  suppressWarnings({
    summary.NOTannotateExon$top5OfftargetTotalScore <- as.numeric(summary.NOTannotateExon$top5OfftargetTotalScore)
    summary.NOTannotateExon$top10OfftargetTotalScore <- as.numeric(summary.NOTannotateExon$top10OfftargetTotalScore)
    offtarget.NOTannotateExon$mismatch.distance2PAM <- offtarget.NOTannotateExon$mismatch.distance2PAM %>% as.numeric() %>% as.character()
  })
  
  expect_equal(summary.NOTannotateExon, summary[names(summary.NOTannotateExon)], tolerance = 0.01)
  expect_equal(offtarget.NOTannotateExon, offtarget[names(offtarget.NOTannotateExon)], tolerance = 0.01)
})