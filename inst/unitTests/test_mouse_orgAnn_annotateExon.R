library(CRISPRseek)
library(org.Mm.eg.db)
library("BSgenome.Mmusculus.UCSC.mm10")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
inputFilePath = system.file("extdata", "RIPK1stop.fa", package = "CRISPRseek")
## org.Dm.egFLYBASE2EG
## org.At.tair.db
## org.At.tairSYMBOL
## org.Ce.egSYMBOL
## org.Mm.egSYMBOL
## org.Rn.egSYMBOL
## org.Dr.egSYMBOL
## org.Hs.egSYMBOL

test.gRNAPlusPAM <- TRUE

outputDir <- getwd();
REpatternFile <- system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek")

summary.NOorgAnn <- read.table(
    system.file("extdata/testMouse/NOorgAnn",
    "Summary.xls", package = "CRISPRseek"), sep = "\t", header = TRUE,
    stringsAsFactors = FALSE)
offtarget.NOorgAnn  <- read.table(
    system.file("extdata/testMouse/NOorgAnn",
    "OfftargetAnalysis.xls", package = "CRISPRseek"), sep = "\t",
    header = TRUE, stringsAsFactors = FALSE)

summary.WithorgAnn <- read.table(
	system.file("extdata/testMouse/WithorgAnn",
	"Summary.xls", package = "CRISPRseek"), sep = "\t", header = TRUE,
	stringsAsFactors = FALSE)
offtarget.WithorgAnn  <- read.table(
	system.file("extdata/testMouse/WithorgAnn",
	"OfftargetAnalysis.xls", package = "CRISPRseek"), sep = "\t",
	header = TRUE, stringsAsFactors = FALSE)

summary.NOTannotateExon <- read.table(
	system.file("extdata/testMouse/NOTannotateExon",
	"Summary.xls", package = "CRISPRseek"), sep = "\t", header = TRUE,
	stringsAsFactors = FALSE)
offtarget.NOTannotateExon  <- read.table(
	system.file("extdata/testMouse/NOTannotateExon",
	"OfftargetAnalysis.xls", package = "CRISPRseek"), sep = "\t",
	header = TRUE, stringsAsFactors = FALSE)

summary.WithorgAnn <- summary.WithorgAnn[,1:21]
summary.NOorgAnn <- summary.NOorgAnn [,1:21]
summary.NOTannotateExon <- summary.NOTannotateExon[,1:21]

if (!test.gRNAPlusPAM)
{
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

test_mouse_orgAnn_annotateExon <- function() {
    cat("Testing for mouse with orgAnn...\n")
    offTargetAnalysis(inputFilePath, findgRNAs = TRUE, 
	  findgRNAsWithREcutOnly = FALSE, findPairedgRNAOnly = FALSE,
	  BSgenomeName = Mmusculus, annotateExon=TRUE, 
          chromToSearch = chroms,
          min.score = 0.5, topN = 100, outputDir = outputDir, 
	  overwrite = TRUE, max.mismatch=1, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene,
	  orgAnn = org.Mm.egSYMBOL)
	
    summary <- read.table("Summary.xls", sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE)
    offtarget <- read.table("OfftargetAnalysis.xls", sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE)
   if (!test.gRNAPlusPAM)
   {
       summary <- summary[, -exclude.sum.col]
       offtarget <- offtarget[, -exclude.oft.col]
   }
    if (checkEquals(summary.WithorgAnn, summary, tolerance = 0.01))
        cat("Summary passed test for orgAnn provided\n")
    else
        cat("Summary failed test for orgAnn provided\n")
    if (checkEquals(offtarget.WithorgAnn, offtarget,tolerance = 0.001))
        cat("off target analysis details passed test for orgAnn provided\n")
    else
        cat("off target anlalysis details failed test for orgAnn provided\n")
    cat("Testing for mouse without orgAnn...\n")
    offTargetAnalysis(inputFilePath, findgRNAs = TRUE, 
	  findgRNAsWithREcutOnly = FALSE, findPairedgRNAOnly = FALSE,
	  BSgenomeName = Mmusculus, annotateExon=TRUE, 
          chromToSearch = chroms,
          min.score = 0.5, topN = 100, outputDir = outputDir, 
	  overwrite = TRUE, max.mismatch=1, 
	  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene	)
	
    summary <- read.table("Summary.xls", sep = "\t", header = TRUE, 
						  stringsAsFactors = FALSE)
    offtarget <- read.table("OfftargetAnalysis.xls", sep = "\t", header = TRUE, 
							stringsAsFactors = FALSE)
    if (!test.gRNAPlusPAM)
   {
       summary <- summary[, -exclude.sum.col]
       offtarget <- offtarget[, -exclude.oft.col]
   }
    if (checkEquals(summary.NOorgAnn, summary, tolerance = 0.01))
	cat("Summary passed test for orgAnn not provided\n")
    else
	cat("Summary failed test for orgAnn not provided\n")
    if (checkEquals(offtarget.NOorgAnn, offtarget,tolerance = 0.001))
	cat("off target analysis details passed test for orgAnn not provided\n")
    else
	cat("off target anlalysis details failed test for orgAnn not provided\n")
    cat("Testing for mouse without orgAnn...\n")
    offTargetAnalysis(inputFilePath, findgRNAs = TRUE, 
  	findgRNAsWithREcutOnly = FALSE, findPairedgRNAOnly = FALSE,
        BSgenomeName = Mmusculus, annotateExon=FALSE, 
        chromToSearch = chroms,
         min.score = 0.5, topN = 100, outputDir = outputDir, 
	  overwrite = TRUE, max.mismatch=1)
	
    summary <- read.table("Summary.xls", sep = "\t", header = TRUE, 
	stringsAsFactors = FALSE)
    offtarget <- read.table("OfftargetAnalysis.xls", sep = "\t", header = TRUE, 
	stringsAsFactors = FALSE)
    if (!test.gRNAPlusPAM)
   {
       summary <- summary[, -exclude.sum.col]
       offtarget <- offtarget[, -exclude.oft.col]
   }
    if (checkEquals(summary.NOTannotateExon, summary, tolerance = 0.01))
	cat("Summary passed test for NOT annotateExon\n")
    else
	cat("Summary failed test for NOT annotateExon\n")
    if (checkEquals(offtarget.NOTannotateExon, offtarget,tolerance = 0.001))
	cat("off target analysis details passed test for NOT annotateExon\n")
    else
	cat("off target anlalysis details failed test for NOT annotateExon\n")
}
