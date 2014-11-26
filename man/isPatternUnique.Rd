\name{isPatternUnique}
\alias{isPatternUnique}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Output whether the input patterns occurs only once in the sequence
}
\description{
Input a sequence and a list of patterns and determine if the patterns occurs only once in the sequence. Used for determining whether a RE site in gRNA also occurs in the flanking region.
}
\usage{
isPatternUnique(seq, patterns)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seq}{
flanking sequence of a gRNA
}
  \item{patterns}{
patterns as DNAStringSet, such as a list of RE sites 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 returns a character vectors containing the uniqueness of each pattern/RE site
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
}
\examples{
    seq <- "TGGATTGTATAATCAGCATGGATTTGGAAC"
    patterns <- DNAStringSet(c("TGG", "TGGA", "TGGATA", "TTGGAAC", ""))
    isPatternUnique(seq, patterns)
    isPatternUnique(seq)
    isPatternUnique(patterns)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
