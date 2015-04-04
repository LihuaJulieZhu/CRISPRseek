\name{foldgRNAs}
\alias{foldgRNAs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fold gRNAs with the gRNA backbone constant region
}
\description{
Fold gRNAs with the gRNA backbone constant region and output minimum free energy and the folded structure in bracket notation using GeneRfold package 
}
\usage{
foldgRNAs(gRNAs.withoutPAM, 
    gRNA.backbone="GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU",
    temperature = 37)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gRNAs.withoutPAM}{
gRNAs as character, without PAM sequence.  
}
  \item{gRNA.backbone}{
gRNA backbone constant region sequence. Default to the
sequence in Sp gRNA backbone.
}
  \item{temperature}{
temperature in celsius. Default to 37 celsius. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 a data frame that contains a set of gRNAs 
 annotated with bracket.notation (folded sgRNA in bracket notation),
 mfe.sgRNA (minimum free energy of sgRNA:gRNA plus backbone),
 mfe.diff (mfe.sgRNA-mfe.backbone), mfe.backbone (minimum free energy
 of the gRNA backbone by itself).
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
    gRNAs.withoutPAM <- c("AACCTTGGGGTTACTGAAAG", "ATCCTGGAGCTTAGTCATAG")
    mfe <- foldgRNAs(gRNAs.withoutPAM)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
