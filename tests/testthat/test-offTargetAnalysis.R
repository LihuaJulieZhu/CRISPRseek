test_that("test_offTargetAnalysis", {
    outputDir <- getwd();
    chroms <- c("chr22", "chrX", "chr15", "chr18", "chr20", 
                "chr5", "chr17", "chr19", "chr16", "chr10", "chr9", "chr1", "chr8", "chr12", "chrY", "chr2", "chr14", "chr11", "chr3", "chr7", "chr13", "chr4", "chr6")


    test.gRNAPlusPAM <-FALSE

    inputFilePath <- system.file("extdata", "inputseq.fa", package = "CRISPRseek")
    REpatternFile <- system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek")
    gRNAFilePath <- system.file("extdata", "testHsap_GATA1_ex2_gRNA1.fa",
     package = "CRISPRseek")
    REcutDetails.pairedSearchRE <- 
        read.table(system.file("extdata/pairedSearch/withRE", 
        "REcutDetails.xls", package = "CRISPRseek"), sep = "\t", 
        header = TRUE, stringsAsFactors = FALSE)
    pairedgRNAs.pairedSearchRE <- read.table(
        system.file("extdata/pairedSearch/withRE",
    "pairedgRNAs.xls", package = "CRISPRseek"), sep = "\t",
        header = TRUE, stringsAsFactors = FALSE)
    summary.pairedSearchRE <- read.table(
        system.file("extdata/pairedSearch/withRE",
        "Summary.xls", package = "CRISPRseek"), sep = "\t", header = TRUE,
        stringsAsFactors = FALSE)
    offtarget.pairedSearchRE  <- read.table(
        system.file("extdata/pairedSearch/withRE",
        "OfftargetAnalysis.xls", package = "CRISPRseek"), sep = "\t",
        header = TRUE, stringsAsFactors = FALSE)

    REcutDetails.pairedSearchNoRE <- read.table(
        system.file("extdata/pairedSearch/withoutRE",
        "REcutDetails.xls", package = "CRISPRseek"), sep = "\t",
     header = TRUE, stringsAsFactors = FALSE)
    pairedgRNAs.pairedSearchNoRE <- read.table(
        system.file("extdata/pairedSearch/withoutRE",
        "pairedgRNAs.xls",package = "CRISPRseek"), sep = "\t",
        header = TRUE, stringsAsFactors = FALSE)
    summary.pairedSearchNoRE <- read.table(
        system.file("extdata/pairedSearch/withoutRE",
        "Summary.xls", package = "CRISPRseek"), sep = "\t",
        header = TRUE, stringsAsFactors = FALSE)
    offtarget.pairedSearchNoRE  <- read.table(
        system.file("extdata/pairedSearch/withoutRE",
        "OfftargetAnalysis.xls", package = "CRISPRseek"), sep = "\t",
        header = TRUE, stringsAsFactors = FALSE)

    REcutDetails.unPairedSearchRE <- read.table(
        system.file("extdata/unPairedSearch/withRE",
        "REcutDetails.xls", package = "CRISPRseek"), sep = "\t",
        header = TRUE, stringsAsFactors = FALSE)
    pairedgRNAs.unPairedSearchRE <- read.table(
        system.file("extdata/unPairedSearch/withRE",
        "pairedgRNAs.xls", package = "CRISPRseek"), sep = "\t",
        header = TRUE, stringsAsFactors = FALSE)
    summary.unPairedSearchRE <- read.table(
        system.file("extdata/unPairedSearch/withRE",
        "Summary.xls", package = "CRISPRseek"), sep = "\t", header = TRUE,
        stringsAsFactors = FALSE)
    offtarget.unPairedSearchRE  <- read.table(
        system.file("extdata/unPairedSearch/withRE",
        "OfftargetAnalysis.xls", package = "CRISPRseek"), sep = "\t", 
        header = TRUE, stringsAsFactors = FALSE)

    REcutDetails.unPairedSearchNoRE <- read.table(
        system.file("extdata/unPairedSearch/withoutRE",
        "REcutDetails.xls", package = "CRISPRseek"), sep = "\t", 
     header = TRUE, stringsAsFactors = FALSE)
    pairedgRNAs.unPairedSearchNoRE <- read.table(
        system.file("extdata/unPairedSearch/withoutRE",
        "pairedgRNAs.xls", package = "CRISPRseek"), sep = "\t",
        header = TRUE, stringsAsFactors = FALSE)
    summary.unPairedSearchNoRE <- read.table(
        system.file("extdata/unPairedSearch/withoutRE",
        "Summary.xls",package = "CRISPRseek"), sep = "\t",
        header = TRUE, stringsAsFactors = FALSE)
    offtarget.unPairedSearchNoRE  <- read.table(
        system.file("extdata/unPairedSearch/withoutRE",
        "OfftargetAnalysis.xls", package = "CRISPRseek"), sep = "\t",
         header = TRUE, stringsAsFactors = FALSE)

    REcutDetails.gRNAProvided <- read.table(
        system.file("extdata/gRNAProvidedSearch",
        "REcutDetails.xls", package = "CRISPRseek"), sep = "\t",
        header = TRUE, stringsAsFactors = FALSE)
    summary.gRNAProvided <- read.table(
        system.file("extdata/gRNAProvidedSearch",
        "Summary.xls", package = "CRISPRseek"), sep = "\t", header = TRUE,
        stringsAsFactors = FALSE)
    offtarget.gRNAProvided  <- read.table(
        system.file("extdata/gRNAProvidedSearch",
        "OfftargetAnalysis.xls", package = "CRISPRseek"), sep = "\t",
        header = TRUE, stringsAsFactors = FALSE)

    summary.gRNAProvided <- summary.gRNAProvided[,1:21]
    summary.unPairedSearchNoRE <- summary.unPairedSearchNoRE[,1:21]
    summary.pairedSearchRE <- summary.pairedSearchRE[,1:22]
    summary.pairedSearchNoRE <- summary.pairedSearchNoRE[,1:22]
    summary.unPairedSearchRE <- summary.unPairedSearchRE[,1:21]

    if (!test.gRNAPlusPAM)
    {
        exclude.sum.col <-  grep("gRNAsPlusPAM", colnames(summary.gRNAProvided))
        exclude.oft.col <- grep("gRNAPlusPAM", colnames(offtarget.gRNAProvided))
        offtarget.gRNAProvided <- offtarget.gRNAProvided[, -exclude.oft.col]
        summary.gRNAProvided <- summary.gRNAProvided[, -exclude.sum.col]
        offtarget.unPairedSearchNoRE <- offtarget.unPairedSearchNoRE[, -exclude.oft.col]
        summary.unPairedSearchNoRE <- summary.unPairedSearchNoRE[, -exclude.sum.col]
        summary.pairedSearchRE <- summary.pairedSearchRE[, -exclude.sum.col]
        offtarget.pairedSearchRE <- offtarget.pairedSearchRE[, -exclude.oft.col]
        offtarget.unPairedSearchRE <- offtarget.unPairedSearchRE[, -exclude.oft.col]
        summary.unPairedSearchRE <- summary.unPairedSearchRE[, -exclude.sum.col]
        summary.pairedSearchNoRE <- summary.pairedSearchNoRE[, -exclude.sum.col]
        offtarget.pairedSearchNoRE <- offtarget.pairedSearchNoRE[, -exclude.oft.col]
    }
    cat("Testing for findgRNAs = FALSE...\n")
    offTargetAnalysis(inputFilePath = gRNAFilePath, findgRNAs = FALSE, 
        findgRNAsWithREcutOnly = FALSE, REpatternFile = REpatternFile, 
        findPairedgRNAOnly = FALSE, BSgenomeName = Hsapiens, 
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
	orgAnn = org.Hs.egSYMBOL, max.mismatch = 3, 
        chromToSearch = c("chrX", "chr11", "chr4", "chr5", "chr6", "chr8", "chr15", "chr3"),
        min.score = 0.5, topN = 100,
        outputDir = outputDir, overwrite = TRUE)
    REcutDetails <- read.table("REcutDetails.xls", sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE)
    summary <- read.table("Summary.xls", sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE)
    offtarget <- read.table("OfftargetAnalysis.xls", sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE)

    if (!test.gRNAPlusPAM)
   {
       summary <- summary[, -exclude.sum.col]
       offtarget <- offtarget[, -exclude.oft.col]
   }

    expect_equal(REcutDetails.gRNAProvided, REcutDetails)
    expect_equal(summary.gRNAProvided, summary, tolerance = 0.01)
    expect_equal(offtarget.gRNAProvided, offtarget,tolerance = 0.001)
    for (isPaired in c(TRUE, FALSE))
    {	
        for (isRE in c(TRUE, FALSE))
        {
            offTargetAnalysis(inputFilePath = inputFilePath, 
                findgRNAsWithREcutOnly = isRE, REpatternFile = REpatternFile,
                findPairedgRNAOnly = isPaired, BSgenomeName = Hsapiens, 
                txdb=TxDb.Hsapiens.UCSC.hg19.knownGene, 
                chromToSearch = chroms,
		orgAnn = org.Hs.egSYMBOL, max.mismatch = 3, 
                min.score = 0.5, topN = 100,
                outputDir= outputDir, overwrite = TRUE)
            REcutDetails <- read.table("REcutDetails.xls", sep = "\t", 
                header = TRUE, stringsAsFactors = FALSE)
            pairedgRNAs <- read.table("pairedgRNAs.xls", sep = "\t", 
                header = TRUE, stringsAsFactors = FALSE)
            summary <- read.table("Summary.xls", sep = "\t", header = TRUE, 
                stringsAsFactors = FALSE)
            offtarget <- read.table("OfftargetAnalysis.xls", sep = "\t", 
                header = TRUE, stringsAsFactors = FALSE)
            if (!test.gRNAPlusPAM)
            {
                summary <- summary[, -exclude.sum.col]
                offtarget <- offtarget[, -exclude.oft.col]
             }

            if (isPaired && isRE)
            {
                cat("Testing for paired with RE sites...\n")
                expect_equal(REcutDetails.pairedSearchRE, REcutDetails )
                expect_equal(pairedgRNAs.pairedSearchRE, pairedgRNAs)
                expect_equal(
                    summary.pairedSearchRE, summary, tolerance = 0.01)
                expect_equal(
                    offtarget.pairedSearchRE, offtarget, tolerance = 0.001)
            }
            else if (isPaired && !isRE)
            {
                cat("Testing for paired with or without RE sites...\n")
                expect_equal(REcutDetails.pairedSearchNoRE, REcutDetails )
                expect_equal(pairedgRNAs.pairedSearchNoRE, pairedgRNAs)
                expect_equal(
                    summary.pairedSearchNoRE, summary, tolerance = 0.01)
                expect_equal(
                    offtarget.pairedSearchNoRE, offtarget, tolerance = 0.001)
            }
            else if (!isPaired && isRE)
            {
                cat("Testing for unPaired with RE sites...\n")
                expect_equal(REcutDetails.unPairedSearchRE, REcutDetails )
                expect_equal(pairedgRNAs.unPairedSearchRE, pairedgRNAs)
                expect_equal(
                    summary.unPairedSearchRE, summary, tolerance = 0.01)
                expect_equal(
                    offtarget.unPairedSearchRE, offtarget, tolerance = 0.001)
            }
            else
            {
                cat("Testing for unPaired with or without RE sites...\n")
                expect_equal(REcutDetails.unPairedSearchNoRE, REcutDetails )
                expect_equal(
                    pairedgRNAs.unPairedSearchNoRE, pairedgRNAs)
                expect_equal(
                    summary.unPairedSearchNoRE, summary, tolerance = 0.01)
                expect_equal(
                    offtarget.unPairedSearchNoRE, offtarget, tolerance = 0.001)
            }
        }
    }
})
