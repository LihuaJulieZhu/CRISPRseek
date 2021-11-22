test_that("test_gRNAprovidedNoPAM", {
    cat("testing gRNAs provided and filter gRNAs with RE for offTargetAnalysis\n")
    inputgRNA <- DNAStringSet("GTAGATGAGGGAGCAGGCGTTGG")
    names(inputgRNA) <- "rs362331T"
    chroms <- "chr4"
    gRNAprovidePAM <- offTargetAnalysis(inputgRNA, 
       findgRNAsWithREcutOnly = FALSE, gRNAoutputName = "testgRNAprovideNoPAM",
       findPairedgRNAOnly = FALSE, findgRNAs = FALSE,
       annotatePaired = FALSE,
       BSgenomeName = Hsapiens, chromToSearch = chroms,
       max.mismatch = 1, annotateExon = FALSE,
       outputDir = getwd(), overwrite = TRUE)

    gRNAprovidePAMre <- offTargetAnalysis(inputgRNA, 
        findgRNAsWithREcutOnly = TRUE, gRNAoutputName = "testgRNAprovideNoPAM",
        findPairedgRNAOnly = FALSE, findgRNAs = FALSE,
        annotatePaired = FALSE,
        BSgenomeName = Hsapiens, chromToSearch = chroms,
        max.mismatch = 1, annotateExon = FALSE,
        outputDir = getwd(), overwrite = TRUE)

    inputgRNA <- DNAStringSet("GTAGATGAGGGAGCAGGCGT")
    names(inputgRNA) <- "rs362331T"

    gRNAprovideNoPAM <- offTargetAnalysis(inputgRNA, 
       findgRNAsWithREcutOnly = FALSE, gRNAoutputName = "testgRNAprovideNoPAM",
       findPairedgRNAOnly = FALSE, findgRNAs = FALSE,
       annotatePaired = FALSE,
       BSgenomeName = Hsapiens, chromToSearch = chroms,
       max.mismatch = 1, annotateExon = FALSE,
       outputDir = getwd(), overwrite = TRUE)

    gRNAprovideNoPAMre <- offTargetAnalysis(inputgRNA,
        findgRNAsWithREcutOnly = TRUE, gRNAoutputName = "testgRNAprovideNoPAM",
        findPairedgRNAOnly = FALSE, findgRNAs = FALSE,
        annotatePaired = FALSE,
        BSgenomeName = Hsapiens, chromToSearch = chroms,
        max.mismatch = 1, annotateExon = FALSE,
        outputDir = getwd(), overwrite = TRUE)

    expect_equal(gRNAprovidePAMre$offtarget[, -2], gRNAprovideNoPAMre$offtarget[, 
        -2]) 
    cat("testing gRNAs provided without filtering gRNAs with RE for offTargetAnalysis\n")
    expect_equal(gRNAprovidePAM$offtarget[, -2], gRNAprovideNoPAM$offtarget[, 
        -2])
})
