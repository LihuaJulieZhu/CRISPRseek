test_that("test_gRNAprovidedCompare2Seqs", {
    inputFile2Path <- system.file("extdata", "rs362331C.fa",
        package = "CRISPRseek")
    seqs2 <- readRDS(system.file("extdata", "seqs2.RDS",
        package = "CRISPRseek"))

    inputgRNA <- DNAStringSet("GTAGATGAGGGAGCAGGCGT")

    gRNAs1Provided <- compare2Sequences(inputgRNA, inputFile2Path,
        outputDir = getwd(), findgRNAs = c(FALSE, FALSE), searchDirection = "1to2",
         overwrite = TRUE)
    gRNAs2Provided <- compare2Sequences( inputFile2Path, inputgRNA,
        outputDir = getwd(), findgRNAs = c(FALSE, FALSE), searchDirection = "2to1",
         overwrite = TRUE)

    for (i in 2:11) {
       expect_equal(as.character(gRNAs1Provided[1, i]), 
            as.character(seqs2[1, i]))
    }
    expect_equal(gRNAs1Provided[, -c(3, 4, 7, 8, 14)], gRNAs2Provided[, 
        -c(3, 4, 7, 8, 14)])
    expect_equal(as.character(gRNAs1Provided[, 3]), as.character(gRNAs2Provided[, 
        4])) 
    expect_equal(as.character(gRNAs1Provided[, 4]), as.character(gRNAs2Provided[, 
        3]))
    expect_equal(as.character(gRNAs1Provided[, 7]), as.character(gRNAs2Provided[, 
        8]))
    expect_equal(as.character(gRNAs1Provided[, 8]), as.character(gRNAs2Provided[, 
        7]))
    expect_equal(abs(as.numeric(as.character(gRNAs1Provided[, 
        14]))), abs(as.numeric(as.character(gRNAs2Provided[, 
        14]))))
})

