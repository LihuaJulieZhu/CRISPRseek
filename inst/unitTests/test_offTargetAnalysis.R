library(CRISPRseek)
library("BSgenome.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
outputDir <- getwd();
inputFilePath <- system.file("extdata", "inputseq.fa", package = "CRISPRseek")
REpatternFile <- system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek")
spacerFilePath <- system.file("extdata", "testHsap_GATA1_ex2_gRNA1.fa",
     package = "CRISPRseek")
REcutDetails.pairedSearchRE <- read.table(system.file("extdata/pairedSearch/withRE", 
    "REcutDetails.xls", package = "CRISPRseek"), sep = "\t", 
    header = TRUE, stringsAsFactors = FALSE)
pairedSpacers.pairedSearchRE <- read.table(
    system.file("extdata/pairedSearch/withRE",
    "pairedSpacers.xls", package = "CRISPRseek"), sep = "\t",
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
pairedSpacers.pairedSearchNoRE <- read.table(
    system.file("extdata/pairedSearch/withoutRE",
    "pairedSpacers.xls",package = "CRISPRseek"), sep = "\t",
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
pairedSpacers.unPairedSearchRE <- read.table(
    system.file("extdata/unPairedSearch/withRE",
    "pairedSpacers.xls", package = "CRISPRseek"), sep = "\t",
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
pairedSpacers.unPairedSearchNoRE <- read.table(
    system.file("extdata/unPairedSearch/withoutRE",
    "pairedSpacers.xls", package = "CRISPRseek"), sep = "\t",
     header = TRUE, stringsAsFactors = FALSE)
summary.unPairedSearchNoRE <- read.table(
    system.file("extdata/unPairedSearch/withoutRE",
    "Summary.xls",package = "CRISPRseek"), sep = "\t",
    header = TRUE, stringsAsFactors = FALSE)
offtarget.unPairedSearchNoRE  <- read.table(
    system.file("extdata/unPairedSearch/withoutRE",
    "OfftargetAnalysis.xls", package = "CRISPRseek"), sep = "\t",
    header = TRUE, stringsAsFactors = FALSE)

REcutDetails.spacerProvided <- read.table(
    system.file("extdata/spacerProvidedSearch",
    "REcutDetails.xls", package = "CRISPRseek"), sep = "\t",
    header = TRUE, stringsAsFactors = FALSE)
summary.spacerProvided <- read.table(
    system.file("extdata/spacerProvidedSearch",
    "Summary.xls", package = "CRISPRseek"), sep = "\t", header = TRUE,
    stringsAsFactors = FALSE)
offtarget.spacerProvided  <- read.table(
    system.file("extdata/spacerProvidedSearch",
    "OfftargetAnalysis.xls", package = "CRISPRseek"), sep = "\t",
    header = TRUE, stringsAsFactors = FALSE)

test_offTargetAnalysis <- function() {
    cat("Testing for findSpacers = FALSE...\n")
    offTargetAnalysis(inputFilePath = spacerFilePath, findSpacers = FALSE, 
        findSpacersWithREcutOnly = FALSE, REpatternFile = REpatternFile, 
        findPairedSpacerOnly = FALSE, BSgenomeName = Hsapiens, 
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, max.mismatch = 3, 
        outputDir = outputDir, overwrite = TRUE)
    REcutDetails <- read.table("REcutDetails.xls", sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE)
    summary <- read.table("Summary.xls", sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE)
    offtarget <- read.table("OfftargetAnalysis.xls", sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE)
    if (checkEquals(REcutDetails.spacerProvided, REcutDetails ))
        cat("REcutDetails passed test for spacer provided\n")
    else
        cat("REcutDetails failed for spacer provided\n")
    if (checkEquals(summary.spacerProvided, summary, tolerance = 0.01))
        cat("Summary passed test for spacer provided\n")
    else
        cat("Summary failed for spacer provided\n")
    if (checkEquals(offtarget.spacerProvided, offtarget,tolerance = 0.001))
        cat("off target analysis details passed test for spacer provided\n")
    else
        cat("off target anlalysis details failed for spacer provided\n")
    for (isPaired in c(TRUE, FALSE))
    {	
        for (isRE in c(TRUE, FALSE))
        {
            offTargetAnalysis(inputFilePath = inputFilePath, 
                findSpacersWithREcutOnly = isRE, REpatternFile = REpatternFile,
                findPairedSpacerOnly = isPaired, BSgenomeName = Hsapiens, 
                txdb=TxDb.Hsapiens.UCSC.hg19.knownGene, max.mismatch = 3, 
                outputDir= outputDir, overwrite = TRUE)
            REcutDetails <- read.table("REcutDetails.xls", sep = "\t", 
                header = TRUE, stringsAsFactors = FALSE)
            pairedSpacers <- read.table("pairedSpacers.xls", sep = "\t", 
                header = TRUE, stringsAsFactors = FALSE)
            summary <- read.table("Summary.xls", sep = "\t", header = TRUE, 
                stringsAsFactors = FALSE)
            offtarget <- read.table("OfftargetAnalysis.xls", sep = "\t", 
                header = TRUE, stringsAsFactors = FALSE)
            if (isPaired && isRE)
            {
                cat("Testing for paired with RE sites...\n")
                if (checkEquals(REcutDetails.pairedSearchRE, REcutDetails ))
                    cat("REcutDetails passed test for paired with RE sites\n")
                else
                    cat("REcutDetails failed for paired with RE sites\n")
                if (checkEquals(pairedSpacers.pairedSearchRE, pairedSpacers))
                    cat("pairedSpacers passed test for paired with RE sites\n")
                else
                    cat("pairedSpacers failed for paired with RE sites\n")	
                if (checkEquals(
                    summary.pairedSearchRE, summary, tolerance = 0.01))
                    cat("Summary passed test for paired with RE sites\n")
                else
                    cat("Summary failed for paired with RE sites\n")
                if (checkEquals(
                    offtarget.pairedSearchRE, offtarget, tolerance = 0.001))
                    cat("off target details passed test for paired 
                        with RE sites\n")
                else
                    cat("off target details failed for paired with RE sites\n")
            }
            else if (isPaired && !isRE)
            {
                cat("Testing for paired with or without RE sites...\n")
                if (checkEquals(REcutDetails.pairedSearchNoRE, REcutDetails ))
                    cat("REcutDetails passed test for paired with or without 
                        RE sites\n")
                else
                    cat("REcutDetails failed for paired with or without RE 
                        sites\n")
                if (checkEquals(pairedSpacers.pairedSearchNoRE, pairedSpacers))
                    cat("pairedSpacers passed test for paired with or without 
                        RE sites\n")
                else
                    cat("pairedSpacers failed for paired with or without RE 
                        sites\n")
                if (checkEquals(
                    summary.pairedSearchNoRE, summary, tolerance = 0.01))
                    cat("Summary passed test for paired with or without RE 
                        sites\n")
                else
                    cat("Summary failed\n")
                if (checkEquals(
                    offtarget.pairedSearchNoRE, offtarget, tolerance = 0.001))
                    cat("off target details passed test for paired with or 
                        without RE sites\n")
                else
                    cat("off target details failed for paired with or without 
                        RE sites\n")
            }
            else if (!isPaired && isRE)
            {
                cat("Testing for unPaired with RE sites...\n")
                if (checkEquals(REcutDetails.unPairedSearchRE, REcutDetails ))
                    cat("REcutDetails passed test for unPaired with RE sites\n")
                else
                    cat("REcutDetails failed for unPaired with RE sites\n")
                if (checkEquals(pairedSpacers.unPairedSearchRE, pairedSpacers))
                    cat("pairedSpacers passed test for unPaired with RE sites\n")
                else
                    cat("pairedSpacers failed for unPaired with RE sites\n")
                if (checkEquals(
                    summary.unPairedSearchRE, summary, tolerance = 0.01))
                    cat("Summary passed test for unPaired with RE sites\n")
                else
                    cat("Summary failed for unPaired with RE sites\n")
                if (checkEquals(
                    offtarget.unPairedSearchRE, offtarget, tolerance = 0.001))
                    cat("off target details passed test for unPaired with RE 
                        sites\n")
                else
                    cat("off target details failed for unPaired with RE 
                        sites\n")
            }
            else
            {
                cat("Testing for unPaired with or without RE sites...\n")
                if (checkEquals(REcutDetails.unPairedSearchNoRE, REcutDetails ))
                    cat("REcutDetails passed test unPaired with or without RE 
                        sites\n")
                else
                    cat("REcutDetails failed unPaired with or without RE 
                        sites\n")
                if (checkEquals(
                    pairedSpacers.unPairedSearchNoRE, pairedSpacers))
                    cat("pairedSpacers passed test unPaired with or without RE 
                        sites\n")
                else
                    cat("pairedSpacers failed unPaired with or without RE 
                        sites\n")
                if (checkEquals(
                    summary.unPairedSearchNoRE, summary, tolerance = 0.01))
                    cat("Summary passed test unPaired with or without RE 
                        sites\n")
                else
                    cat("Summary failed unPaired with or without RE sites\n")
                if (checkEquals(
                    offtarget.unPairedSearchNoRE, offtarget, tolerance = 0.001))
                    cat("off target details passed test unPaired with or 
                        without RE sites\n")
                else
                    cat("off target details failed unPaired with or without RE 
                        sites\n")
            }
        }
    }
}
