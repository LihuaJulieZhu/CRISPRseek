\name{annotateOffTargets}
\alias{annotateOffTargets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
annotate off targets
}
\description{
 annotate off targets to indicate whether it is inside an exon or intron, and the gene id if inside the gene.
}
\usage{
annotateOffTargets(scores, txdb, orgAnn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scores}{
a data frame output from getOfftargetScore or filterOfftarget. It contains strand (strand of the 
off target, + for plus and - for minus strand), chrom (chromosome of the off 
target), chromStart (start position of the off target), chromEnd (end position 
of the off target),name (gRNA name),gRNAPlusPAM (gRNA sequence with PAM 
sequence concatenated), OffTargetSequence (the genomic sequence of the off 
target), n.mismatch (number of mismatches between the off target and the 
gRNA), forViewInUCSC (string for viewing in UCSC genome browser, e.g., 
chr14:31665685-31665707), score (score of the off target), 
mismatch.distance2PAM (a comma separated distances of all mismatches to 
PAM, e.g., 14,11 means one mismatch is 14 bp away from PAM and the other 
mismatch is 11 bp away from PAM), alignment (alignment between gRNA and off 
target, e.g., ......G..C.......... means that this off target aligns with 
gRNA except that G and C are mismatches),NGG (this off target contains 
canonical PAM or not, 1 for yes and 0 for no)  mean.neighbor.distance.mismatch 
(mean distance between neighboring mismatches)
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
 \item{orgAnn}{
organism annotation mapping such as org.Hs.egSYMBOL in org.Hs.eg.db package for human
}
 }
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 a data frame with off target annotation
}
\references{
Lihua Julie Zhu, Benjamin R. Holmes, Neil Aronin and Michael Brodsky. CRISPRseek: a Bioconductor package to identify target-specific guide RNAs for CRISPR-Cas9 genome-editing systems. Plos One Sept 23rd 2014
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
    library(CRISPRseek)
    #library("BSgenome.Hsapiens.UCSC.hg19")
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
    hitsFile <-  system.file("extdata", "hits.txt", package="CRISPRseek")
    hits <- read.table(hitsFile, sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE)
    featureVectors <- buildFeatureVectorForScoring(hits)
    scores <- getOfftargetScore(featureVectors)
    outputDir <- getwd() 
    results <- annotateOffTargets(scores, 
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
         orgAnn = org.Hs.egSYMBOL)
    results
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
