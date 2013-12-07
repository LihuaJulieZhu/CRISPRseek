\name{filterSpacers}
\alias{filterSpacers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter spacers
}
\description{
Filter spacers containing restriction enzyme cut site
}
\usage{
filterSpacers(all.spacers, pairOutputFile = "", REpatternFile, format = "fasta", 
    minREpatternSize = 6, overlap.spacer.positions = c(17, 18))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{all.spacers}{
Spacers as DNAStringSet, such as the output from findSpacers
}
  \item{pairOutputFile}{
File path with paired spacers
}
  \item{REpatternFile}{
File path containing restriction enzyme cut patters
}
  \item{format}{
Format of the REpatternFile, default as fasta
}
  \item{minREpatternSize}{
Minimum restriction enzyme recognition pattern length required for the 
enzyme pattern to be searched for, default 6
}
  \item{overlap.spacer.positions}{
The required overlap positions of spacer and restriction enzyme cut site, 
default 17 and 18
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{spacers.withRE }{spacers as DNAStringSet that passed the filter criteria}
 \item{spacerREcutDetails}{a data frame that contains a set of spacers 
 annotated with restriction enzyme cut details}
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
    all.spacers <- findSpacers(
        inputFilePath = system.file("extdata", "inputseq.fa", 
        package = "CRISPRseek"),
        pairOutputFile = "testpairedSpacers.xls",
        findPairedSpacerOnly = TRUE)

    spacers.RE <- filterSpacers(all.spacers = all.spacers,
        pairOutputFile = "testpairedSpacers.xls",
        REpatternFile = system.file("extdata", "NEBenzymes.fa", 
        package = "CRISPRseek"))

    spacers  <- spacers.RE$spacers.withRE
    restriction.enzyme.cut.sites <- spacers.RE$spacerREcutDetails
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
