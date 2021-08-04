library(CRISPRseek)
library("BSgenome.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library("BSgenome.Mmusculus.UCSC.mm10")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

outputDir <- getwd();

inputFile1Path <- system.file("extdata", "rs362331T.fa",
    package = "CRISPRseek")
inputFile2Path <- system.file("extdata", "rs362331C.fa",
    package = "CRISPRseek")

inputFilePath <- system.file("extdata", "inputseq.fa", package = "CRISPRseek")
REpatternFile <- system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek")
gRNAFilePath <- system.file("extdata", "testHsap_GATA1_ex2_gRNA1.fa",
     package = "CRISPRseek")
outputPairedSearchRE <- "~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata/pairedSearch/withRE"
	
outputPairedSearchNoRE <- "~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata/pairedSearch/withoutRE"
   
outputUnPairedSearchRE <- "~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata/unPairedSearch/withRE"
  
outputUnPairedSearchNoRE <- "~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata/unPairedSearch/withoutRE"

outputgRNAProvided <- "~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata/gRNAProvidedSearch"

### recreate data for notarget found
inputFilePath <- system.file("extdata/testNoTarget/",
    "negativeCntl-gRNA.fa", package = "CRISPRseek")

chroms <- c("chr19", "chr18", "chr13",  "chr12", "chr11")

offTargetAnalysis(inputFilePath,
                   scoring.method = "CFDscore",
                   min.score = 0.001,
                   format = "fasta",
                   findgRNAs = FALSE,
                   findgRNAsWithREcutOnly = FALSE,
                   findPairedgRNAOnly = FALSE,
                   annotatePaired = FALSE,
                   BSgenomeName = Mmusculus,
                   chromToSearch= chroms,
                   max.mismatch = 3,
                   annotateExon = FALSE,
                   fetchSequence = FALSE,
       outputDir = "~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata/testNoTarget",
                   overwrite = TRUE)

## recreate data for compare2sequences
inputFile1Path <- system.file("extdata", "rs362331T.fa",
    package = "CRISPRseek")
inputFile2Path <- system.file("extdata", "rs362331C.fa",
    package = "CRISPRseek")


seqs2CFD.new <- compare2Sequences(inputFile1Path, inputFile2Path,
        outputDir = getwd(),
        overwrite = TRUE,scoring.method = "CFDscore")

seqs2.new <- compare2Sequences(inputFile1Path, inputFile2Path,
        outputDir = getwd(),
        overwrite = TRUE)

seqs2.5prime.new <- compare2Sequences(inputFile1Path, inputFile2Path,
        outputDir = getwd(), PAM.location = "5prime", PAM = "GCT", PAM.pattern = "^NCN",
        overwrite = TRUE, subPAM.position = c(1,2))

saveRDS(seqs2.5prime.new, file ="~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata/seqs2.5prime.RDS")
saveRDS(seqs2.new, file = "~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata/seqs2.RDS")
saveRDS(seqs2CFD.new, file = "~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata/seqs2CFD.RDS")

offTargetAnalysis(inputFilePath = gRNAFilePath, findgRNAs = FALSE,
        findgRNAsWithREcutOnly = FALSE, REpatternFile = REpatternFile,
        findPairedgRNAOnly = FALSE, BSgenomeName = Hsapiens,
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
        orgAnn = org.Hs.egSYMBOL, max.mismatch = 3,
        chromToSearch = c("chrX", "chr11", "chr4", "chr5", "chr6", "chr8", "chr15", "chr3"),
        min.score = 0.5, topN = 100,
        outputDir = outputgRNAProvided, overwrite = TRUE)

   
library("BSgenome.Dmelanogaster.UCSC.dm3")
library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
library("org.Dm.eg.db")
gRNAFilePath <- system.file("extdata", "tobi2.fa", package = "CRISPRseek")
outputNoOffTarget <- "~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata/gRNAProvidedSearchNOoffTarget"
outputNoOffTargetNoAnn <- "~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata/gRNAProvidedSearchNOoffTargetNoAnn"

cat("Creating output for findgRNAs = FALSE and no off target found and annotateExon FALSE...\n")
offTargetAnalysis(inputFilePath = gRNAFilePath, findgRNAs = FALSE,
        findgRNAsWithREcutOnly = FALSE, REpatternFile = REpatternFile,
        findPairedgRNAOnly = FALSE, BSgenomeName = Dmelanogaster,
        txdb = TxDb.Dmelanogaster.UCSC.dm3.ensGene, annotateExon = FALSE,
        orgAnn = org.Dm.egFLYBASE2EG, max.mismatch = 3,
        min.score = 0.5, topN = 100,
        outputDir = outputNoOffTargetNoAnn, overwrite = TRUE)

cat("Creating output for findgRNAs = FALSE and no off target found, annotateExon TRUE...\n")
offTargetAnalysis(inputFilePath = gRNAFilePath, findgRNAs = FALSE,
        findgRNAsWithREcutOnly = FALSE, REpatternFile = REpatternFile,
        findPairedgRNAOnly = FALSE, BSgenomeName = Dmelanogaster,
        txdb = TxDb.Dmelanogaster.UCSC.dm3.ensGene,
        orgAnn = org.Dm.egFLYBASE2EG, max.mismatch = 3,
        min.score = 0.5, topN = 100,
        outputDir = outputNoOffTarget, overwrite = TRUE)

chroms <- c("chr22", "chrX", "chr15", "chr18", "chr20", "chr5", "chr17", "chr19", "chr16", "chr10", "chr9", "chr1", "chr8", "chr12", "chrY", "chr2", "chr14", "chr11", "chr3", "chr7", "chr13", "chr4", "chr6")

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
                chromToSearch = chroms,                
                min.score = 0.5, topN = 100,
                outputDir= outputDir, overwrite = TRUE)
		}
    }


inputFilePath = system.file("extdata", "RIPK1stop.fa", package = "CRISPRseek")
#inputFilePath = "~/CRISPRseekDemo/RIPK1stop.fa"
outputDir = "~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata/testMouse/WithorgAnn"
results <- offTargetAnalysis(inputFilePath, findgRNAs = TRUE, 
findgRNAsWithREcutOnly = FALSE, findPairedgRNAOnly = FALSE,
BSgenomeName = Mmusculus, annotateExon=TRUE, outputDir = outputDir, 
overwrite = TRUE, max.mismatch=1, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene,
min.score = 0.5, topN = 100,
chromToSearch = chroms,
orgAnn = org.Mm.egSYMBOL)

outputDir = "~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata/testMouse/NOorgAnn"
results <- offTargetAnalysis(inputFilePath, findgRNAs = TRUE, 
findgRNAsWithREcutOnly = FALSE, findPairedgRNAOnly = FALSE, 
BSgenomeName = Mmusculus, annotateExon=TRUE, outputDir = outputDir, 
overwrite = TRUE, max.mismatch=1,
min.score = 0.5, topN = 100,
chromToSearch = chroms,
txdb = TxDb.Mmusculus.UCSC.mm10.knownGene)

outputDir = "~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata/testMouse/NOTannotateExon"
results <- offTargetAnalysis(inputFilePath, findgRNAs = TRUE, 
findgRNAsWithREcutOnly = FALSE, findPairedgRNAOnly = FALSE, 
BSgenomeName = Mmusculus, annotateExon=FALSE, outputDir = outputDir, 
min.score = 0.5, topN = 100,
chromToSearch = chroms,
overwrite = TRUE, max.mismatch=1)
