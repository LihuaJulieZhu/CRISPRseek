\name{filterOffTarget}
\alias{filterOffTarget}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
filter off targets and generate reports.
}
\description{
 filter off targets that meet the criteria set by users such as minimum score, 
 topN. In addition, off target was annotated with flank sequence and whether 
 it is inside an exon or not if fetchSequence is set to TRUE and annotateExon 
 is set to TRUE
}
\usage{
filterOffTarget(scores, min.score = 0.5, topN = 100, 
    topN.OfftargetTotalScore = 10,
    annotateExon = TRUE, txdb, outputDir, oneFilePergRNA = FALSE, 
    fetchSequence = TRUE, upstream = 200, downstream = 200, BSgenomeName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scores}{
a data frame output from getOfftargetScore. It contains strand (strand of the 
off target, + for plus and - for minus strand), chrom (chromosome of the off 
target), chromStart (start position of the off target), chromEnd (end position 
of the off target),name (gRNA name),gRNAPlusPAM (gRNA sequence with PAM 
sequence concatenated), OffTargetSequence (the genomic sequence of the off 
target), n.mismatch (number of mismatches between the off target and the 
gRNA), forViewInUCSC (string for viewing in UCSC genome browser, e.g., 
chr14:31665685-31665707), score (score of the off target), 
mismatche.distance2PAM (a comma separated distances of all mismatches to 
PAM, e.g., 14,11 means one mismatch is 14 bp away from PAM and the other 
mismatch is 11 bp away from PAM), alignment (alignment between gRNA and off 
target, e.g., ......G..C.......... means that this off target aligns with 
gRNA except that G and C are mismatches),NGG (this off target contains 
canonical PAM or not, 1 for yes and 0 for no)  mean.neighbor.distance.mismatch 
(mean distance between neighboring mismatches)
}
  \item{min.score}{
minimum score of an off target to included in the final output, default 0.5
}
  \item{topN}{
top N off targets to be included in the final output, default 100
}
  \item{topN.OfftargetTotalScore}{
top N off target used to calculate the total off target score, default 10
}
  \item{annotateExon}{
Choose whether or not to indicate whether the off target is inside an exon or 
not, default TRUE
}
  \item{txdb}{
TxDb object, for creating and using TxDb object, please refer 
to GenomicFeatures package. For a list of existing TxDb object, please
search for annotation package starting with Txdb at 
http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData, 
such as TxDb.Rnorvegicus.UCSC.rn5.refGene for rat, 
TxDb.Mmusculus.UCSC.mm10.knownGene for mouse, 
TxDb.Hsapiens.UCSC.hg19.knownGene for human, 
TxDb.Dmelanogaster.UCSC.dm3.ensGene for Drosophila and 
TxDb.Celegans.UCSC.ce6.ensGene for C.elegans
}
  \item{outputDir}{
the directory where the off target analysis and reports will be written to
}
  \item{oneFilePergRNA}{
write to one file for each gRNA or not, default to FALSE
}
  \item{fetchSequence}{
Fetch flank sequence of off target or not, default TRUE
}
  \item{upstream}{
upstream offset from the off target start, default 200
}
  \item{downstream}{
downstream offset from the off target end, default 200
}
  \item{BSgenomeName}{
BSgenome object. Please refer to available.genomes in BSgenome package. For 
example, BSgenome.Hsapiens.UCSC.hg19 for hg19, BSgenome.Mmusculus.UCSC.mm10 
for mm10, BSgenome.Celegans.UCSC.ce6 for ce6, BSgenome.Rnorvegicus.UCSC.rn5 
for rn5, and BSgenome.Dmelanogaster.UCSC.dm3 for dm3
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{offtargets }{a data frame with off target analysis results}
\item{summary }{a data frame with summary of the off target analysis results}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lihua Julie Zhu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
offTargetAnalysis
}
\examples{
    library("BSgenome.Hsapiens.UCSC.hg19")
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    hitsFile <-  system.file("extdata", "hits.txt", package="CRISPRseek")
    hits <- read.table(hitsFile, sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE)
    featureVectors <- buildFeatureVectorForScoring(hits)
    scores <- getOfftargetScore(featureVectors)
    outputDir <- getwd() 
    results <- filterOffTarget(scores, BSgenomeName = Hsapiens, 
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, outputDir = outputDir, 
        min.score = 0.1, topN = 10, topN.OfftargetTotalScore = 5)
    results$offtargets
    results$summary
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
