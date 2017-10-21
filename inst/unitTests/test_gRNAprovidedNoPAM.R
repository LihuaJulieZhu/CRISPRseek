library("BSgenome.Hsapiens.UCSC.hg19")
library(CRISPRseek)
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

test_gRNAprovidedNoPAM <- function() {
cat("testing gRNAs provided and filter gRNAs with RE for offTargetAnalysis\n")
if(checkEquals(gRNAprovidePAMre$offtarget[,-2] , gRNAprovideNoPAMre$offtarget[,-2]))
  cat("passed test!\n")
else 
  cat("Offtargets of gRNA provided without PAM differs from that of gRNA provided with PAM\n")

cat("testing gRNAs provided without filtering gRNAs with RE for offTargetAnalysis\n")
if(checkEquals(gRNAprovidePAM$offtarget[,-2] , gRNAprovideNoPAM$offtarget[,-2]))
   cat("passed test!\n")
else 
   cat("Offtargets of gRNA provided without PAM differs from that of gRNA provided with PAM\n")
}
