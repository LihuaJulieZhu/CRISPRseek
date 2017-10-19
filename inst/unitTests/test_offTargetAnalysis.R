library(CRISPRseek)
library("BSgenome.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
outputDir <- getwd();

test.gRNAPlusPAM <-FALSE

inputFilePath <- system.file("extdata", "inputseq.fa", package = "CRISPRseek")
REpatternFile <- system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek")
gRNAFilePath <- system.file("extdata", "testHsap_GATA1_ex2_gRNA1.fa",
     package = "CRISPRseek")
REcutDetails.pairedSearchRE <- read.table(system.file("extdata/pairedSearch/withRE", 
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
test_offTargetAnalysis <- function() {
    cat("Testing for findgRNAs = FALSE...\n")
    offTargetAnalysis(inputFilePath = gRNAFilePath, findgRNAs = FALSE, 
        findgRNAsWithREcutOnly = FALSE, REpatternFile = REpatternFile, 
        findPairedgRNAOnly = FALSE, BSgenomeName = Hsapiens, 
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
	orgAnn = org.Hs.egSYMBOL, max.mismatch = 3, 
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

    if (checkEquals(REcutDetails.gRNAProvided, REcutDetails ))
        cat("REcutDetails passed test for gRNA provided\n")
    else
        cat("REcutDetails failed for gRNA provided\n")
    if (checkEquals(summary.gRNAProvided, summary, tolerance = 0.01))
        cat("Summary passed test for gRNA provided\n")
    else
        cat("Summary failed for gRNA provided\n")
    if (checkEquals(offtarget.gRNAProvided, offtarget,tolerance = 0.001))
        cat("off target analysis details passed test for gRNA provided\n")
    else
        cat("off target anlalysis details failed for gRNA provided\n")
    for (isPaired in c(TRUE, FALSE))
    {	
        for (isRE in c(TRUE, FALSE))
        {
            offTargetAnalysis(inputFilePath = inputFilePath, 
                findgRNAsWithREcutOnly = isRE, REpatternFile = REpatternFile,
                findPairedgRNAOnly = isPaired, BSgenomeName = Hsapiens, 
                txdb=TxDb.Hsapiens.UCSC.hg19.knownGene, 
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
                if (checkEquals(REcutDetails.pairedSearchRE, REcutDetails ))
                    cat("REcutDetails passed test for paired with RE sites\n")
                else
                    cat("REcutDetails failed for paired with RE sites\n")
                if (checkEquals(pairedgRNAs.pairedSearchRE, pairedgRNAs))
                    cat("pairedgRNAs passed test for paired with RE sites\n")
                else
                    cat("pairedgRNAs failed for paired with RE sites\n")	
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
                if (checkEquals(pairedgRNAs.pairedSearchNoRE, pairedgRNAs))
                    cat("pairedgRNAs passed test for paired with or without 
                        RE sites\n")
                else
                    cat("pairedgRNAs failed for paired with or without RE 
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
                if (checkEquals(pairedgRNAs.unPairedSearchRE, pairedgRNAs))
                    cat("pairedgRNAs passed test for unPaired with RE sites\n")
                else
                    cat("pairedgRNAs failed for unPaired with RE sites\n")
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
                    pairedgRNAs.unPairedSearchNoRE, pairedgRNAs))
                    cat("pairedgRNAs passed test unPaired with or without RE 
                        sites\n")
                else
                    cat("pairedgRNAs failed unPaired with or without RE 
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
