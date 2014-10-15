library(CRISPRseek)
library("BSgenome.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
outputDir <- getwd();
inputFilePath <- system.file("extdata", "inputseq.fa", package = "CRISPRseek")
REpatternFile <- system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek")
gRNAFilePath <- system.file("extdata", "testHsap_GATA1_ex2_gRNA1.fa",
     package = "CRISPRseek")
outputPairedSearchRE <- "~/Dev/Bioconductor/Trunk/CRISPRseek/inst/extdata/pairedSearch/withRE"
	
outputPairedSearchNoRE <- "~/Dev/Bioconductor/Trunk/CRISPRseek/inst/extdata/pairedSearch/withoutRE"
   
outputUnPairedSearchRE <- "~/Dev/Bioconductor/Trunk/CRISPRseek/inst/extdata/unPairedSearch/withRE"
  
outputUnPairedSearchNoRE <- "~/Dev/Bioconductor/Trunk/CRISPRseek/inst/extdata/unPairedSearch/withoutRE"

outputgRNAProvided <- "~/Dev/Bioconductor/Trunk/CRISPRseek/inst/extdata/gRNAProvidedSearch"

offTargetAnalysis(inputFilePath = gRNAFilePath, findgRNAs = FALSE, 
        findgRNAsWithREcutOnly = FALSE, REpatternFile = REpatternFile, 
        findPairedgRNAOnly = FALSE, BSgenomeName = Hsapiens, 
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
		orgAnn = org.Hs.egSYMBOL, max.mismatch = 3, 
        outputDir = outputgRNAProvided, overwrite = TRUE)
   
library("BSgenome.Dmelanogaster.UCSC.dm3")
library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
library("org.Dm.eg.db")
gRNAFilePath <- system.file("extdata", "tobi2.fa", package = "CRISPRseek")
outputNoOffTarget <- "~/Dev/Bioconductor/Trunk/CRISPRseek/inst/extdata/gRNAProvidedSearchNOoffTarget"
cat("Creating output for findgRNAs = FALSE and no off target found...\n")
offTargetAnalysis(inputFilePath = gRNAFilePath, findgRNAs = FALSE, 
	findgRNAsWithREcutOnly = FALSE, REpatternFile = REpatternFile, 
	findPairedgRNAOnly = FALSE, BSgenomeName = Dmelanogaster, 
	txdb = TxDb.Dmelanogaster.UCSC.dm3.ensGene, 
	orgAnn = org.Dm.egFLYBASE2EG, max.mismatch = 3, 
	outputDir = outputNoOffTarget, overwrite = TRUE)


    for (isPaired in c(TRUE, FALSE))
    {	
        for (isRE in c(TRUE, FALSE))
        {
			if (isRE)
			{
				if (isPaired)
					outputDir <- outputPairedSearchRE
				else
					outputDir <- outputUnPairedSearchRE
			}
			else
			{
				if (isPaired)
					outputDir <- outputPairedSearchNoRE
				else
					outputDir <- outputUnPairedSearchNoRE
			}			
            offTargetAnalysis(inputFilePath = inputFilePath, 
                findgRNAsWithREcutOnly = isRE, REpatternFile = REpatternFile,
                findPairedgRNAOnly = isPaired, BSgenomeName = Hsapiens, 
                txdb=TxDb.Hsapiens.UCSC.hg19.knownGene, 
				orgAnn = org.Hs.egSYMBOL, max.mismatch = 3, 
                outputDir= outputDir, overwrite = TRUE)
		}
    }


library(CRISPRseek)
library(org.Mm.eg.db)
library(CRISPRseek)
library("BSgenome.Mmusculus.UCSC.mm10")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
inputFilePath = system.file("extdata", "RIPK1stop.fa", package = "CRISPRseek")
#inputFilePath = "~/CRISPRseekDemo/RIPK1stop.fa"
outputDir = "~/Dev/Bioconductor/Trunk/CRISPRseek/inst/extdata/testMouse/WithorgAnn"
offTargetAnalysis(inputFilePath, findgRNAs = TRUE, 
findgRNAsWithREcutOnly = FALSE, findPairedgRNAOnly = FALSE,
BSgenomeName = Mmusculus, annotateExon=TRUE, outputDir = outputDir, 
overwrite = TRUE, max.mismatch=1, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene,
orgAnn = org.Mm.egSYMBOL)

outputDir = "~/Dev/Bioconductor/Trunk/CRISPRseek/inst/extdata/testMouse/NOorgAnn"
offTargetAnalysis(inputFilePath, findgRNAs = TRUE, 
findgRNAsWithREcutOnly = FALSE, findPairedgRNAOnly = FALSE, 
BSgenomeName = Mmusculus, annotateExon=TRUE, outputDir = outputDir, 
overwrite = TRUE, max.mismatch=1, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene)

outputDir = "~/Dev/Bioconductor/Trunk/CRISPRseek/inst/extdata/testMouse/NOTannotateExon"
offTargetAnalysis(inputFilePath, findgRNAs = TRUE, 
findgRNAsWithREcutOnly = FALSE, findPairedgRNAOnly = FALSE, 
BSgenomeName = Mmusculus, annotateExon=FALSE, outputDir = outputDir, 
overwrite = TRUE, max.mismatch=1)
