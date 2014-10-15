\name{getOfftargetScore}
\alias{getOfftargetScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate score for each off target
}
\description{
Calculate score for each off target with given feature vectors and weights 
vector
}
\usage{
getOfftargetScore(featureVectors, 
    weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 
    0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{featureVectors}{
a data frame generated from buildFeatureVectorForScoring. It contains 
IsMismatch.posX (Indicator variable indicating whether this position X is 
mismatch or not, 1 means yes and 0 means not, X = 1- gRNA.size) representing 
all positions in the gRNA), strand (strand of the off target, + for plus 
and - for minus strand), chrom (chromosome of the off target), chromStart 
(start position of the off target), chromEnd (end position of the off target),
name (gRNA name),gRNAPlusPAM (gRNA sequence with PAM sequence 
concatenated), OffTargetSequence (the genomic sequence of the off target), 
n.mismatch (number of mismatches between the off target and the gRNA), 
forViewInUCSC (string for viewing in UCSC genome browser, e.g., 
chr14:31665685-31665707), score (score of the off target), 
mismatch.distance2PAM (a comma separated distances of all mismatches to 
PAM, e.g., 14,11 means one mismatch is 14 bp away from PAM and the other 
mismatch is 11 bp away from PAM), alignment (alignment between gRNA and off 
target, e.g., ......G..C.......... means that this off target aligns with 
gRNA except that G and C are mismatches),NGG (this off target contains 
canonical PAM or not, 1 for yes and 0 for no)  mean.neighbor.distance.mismatch 
(mean distance between neighboring mismatches)
}
  \item{weights}{
a numeric vector size of gRNA length, default c(0, 0, 0.014, 0, 0, 0.395, 
0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 
0.685, 0.583) which is used in Hsu et al., 2013 cited in the reference section
}
}
\details{
score is calculated using the weights and algorithm by Hsu et al., 2013 cited 
in the reference section
}
\value{
a data frame containing strand (strand of the match, + for plus and - for minus 
strand), chrom (chromosome of the off target), chromStart (start position of 
the off target), chromEnd (end position of the off target),name (gRNA name), 
gRNAPlusPAM (gRNA sequence with PAM sequence concatenated), 
OffTargetSequence (the genomic sequence of the off target), n.mismatch 
(number of mismatches between the off target and the gRNA), forViewInUCSC 
(string for viewing in UCSC genome browser, e.g., chr14:31665685-31665707), 
score (score of the off target), mismatch.distance2PAM (a comma separated 
distances of all mismatches to PAM, e.g., 14,11 means one mismatch is 14 bp 
away from PAM and the other mismatch is 11 bp away from PAM), alignment 
(alignment between gRNA and off target, e.g., ......G..C.......... means 
that this off target aligns with gRNA except that G and C are mismatches), 
NGG (this off target contains canonical PAM or not, 1 for yes and 0 for no)  
mean.neighbor.distance.mismatch (mean distance between neighboring mismatches)
}
\references{
Patrick D Hsu, David A Scott, Joshua A Weinstein, F Ann Ran, Silvana Konermann, 
Vineeta Agarwala, Yinqing Li, Eli J Fine, Xuebing Wu, Ophir Shalem, Thomas J 
Cradick, Luciano A Marraffini, Gang Bao & Feng Zhang (2013) DNA targeting 
specificity of rNA-guided Cas9 nucleases. Nature Biotechnology 31:827-834
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
    hitsFile <-  system.file("extdata", "hits.txt", 
        package = "CRISPRseek")
    hits <- read.table(hitsFile, sep = "\t", header = TRUE,
        stringsAsFactors = FALSE)
    featureVectors <- buildFeatureVectorForScoring(hits)
    getOfftargetScore(featureVectors)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
