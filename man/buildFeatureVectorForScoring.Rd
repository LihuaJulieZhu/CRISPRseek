\name{buildFeatureVectorForScoring}
\alias{buildFeatureVectorForScoring}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Build feature vectors
}
\description{
Build feature vectors for calculating scores of off targets
}
\usage{
buildFeatureVectorForScoring(hits, gRNA.size = 20, canonical.PAM = "NGG")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hits}{
a data frame generated from searchHits, which contains IsMismatch.posX 
(Indicator variable indicating whether this position X is mismatch or not, 1 
means yes and 0 means not, X = 1- gRNA.size) representing all positions in 
the guide RNA, abbreviated as gRNA),strand (strand of the off target, + for plus and - for minus 
strand),	chrom (chromosome of the off target), chromStart (start position of 
the off target),	chromEnd (end position of the off target),name (gRNA name),	
gRNAPlusPAM (gRNA sequence with PAM sequence concatenated), 
OffTargetSequence (the genomic sequence of the off target), n.mismatch 
(number of mismatches between the off target and the gRNA), forViewInUCSC 
(string for viewing in UCSC genome browser, e.g., chr14:31665685-31665707), 
score (set to 100, and will be calculated in getOfftargetScore)
}
  \item{gRNA.size}{
gRNA size, default 20
}
  \item{canonical.PAM}{
Canonical PAM, default NGG
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    A data frame with hits plus features used for calculating scores and for 
  generating report, including IsMismatch.posX (Indicator variable indicating 
  whether this position X is mismatch or not, 1 means yes and 0 means not, 
  X = 1- gRNA.size) representing all positions in the gRNA), strand 
  (strand of the off target, + for plus and - for minus strand),chrom 
  (chromosome of the off target), chromStart (start position of the off 
  target),chromEnd (end position of the off target),name (gRNA name),	
  gRNAPlusPAM (gRNA sequence with PAM sequence concatenated), 
  OffTargetSequence (the genomic sequence of the off target), n.mismatch 
  (number of mismatches between the off target and the gRNA), forViewInUCSC 
  (string for viewing in UCSC genome browser, e.g., chr14:31665685-31665707), 
  score (score of the off target), mismatche.distance2PAM (a comma separated 
  distances of all mismatches to PAM, e.g., 14,11 means one mismatch is 14 bp 
  away from PAM and the other mismatch is 11 bp away from PAM), alignment 
  (alignment between gRNA and off target, e.g., ......G..C.......... means 
  that this off target aligns with gRNA except that G and C are mismatches), 
  NGG (this off target contains canonical PAM or not, 1 for yes and 0 for no)  
  mean.neighbor.distance.mismatch (mean distance between neighboring mismatches)
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
    hitsFile <-  system.file("extdata", "hits.txt", package = "CRISPRseek")
    hits <- read.table(hitsFile, sep= "\t", header = TRUE,
        stringsAsFactors = FALSE)
    buildFeatureVectorForScoring(hits)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
