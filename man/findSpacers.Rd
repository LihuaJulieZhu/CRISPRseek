\name{findSpacers}
\alias{findSpacers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Find potential spacers
}
\description{
	Find potential spacers for an input file containing sequences in fasta 
    format
}
\usage{
findSpacers(inputFilePath, format = "fasta", PAM = "NGG", PAM.size = 3, 
    findPairedSpacerOnly = FALSE, spacer.size = 20, min.gap = 0, max.gap = 20,
    pairOutputFile, name.prefix = "gRNA")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputFilePath}{
 Sequence input file path that contains sequences to be searched for potential
 spacers
}
  \item{format}{
Format of the input file, fasta and fastq are supported, default fasta
}
  \item{PAM}{
PAM sequence after the spacer, default NGG
}
  \item{PAM.size}{
PAM length, default 3
}
  \item{findPairedSpacerOnly}{
Choose whether to only search for paired spacers in such an orientation that 
the first one is on minus strand called reverse spacer and the second one is 
on plus strand called forward spacer. TRUE or FALSE, default FALSE
}
  \item{spacer.size}{
The size of the spacer, default 20
}
  \item{min.gap}{
Minimum distance between two oppositely oriented spacers to be valid paired 
spacers. Default 0
}
  \item{max.gap}{
Maximum distance between two oppositely oriented spacers to be valid paired 
spacers. Default 20
}
  \item{pairOutputFile}{
The output file for writing paired spacer information to
}
  \item{name.prefix}{
The prefix used when assign name to found spacers, default gRNA, short for
guided RNA.
}
}
\details{
If users already has a fasta file that contains a set of potential spacers, 
then users can call filerSpacers directly although the easiest way is to call 
the one-stop-shopping function OffTargetAnalysis with findSpacers set to FALSE.
}
\value{
DNAStringSet consists of potential spacers that can be input to filterSpacers 
function directly
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lihua Julie Zhu
}
\note{
If the input sequence file contains multiple >300 bp sequences, suggest create
one input file for each sequence and run the OffTargetAnalysis separately.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
offTargetAnalysis
}
\examples{
    findSpacers(inputFilePath = system.file("extdata",
        "inputseq.fa", package = "CRISPRseek"),
        pairOutputFile = "testpairedSpacers.xls",
        findPairedSpacerOnly = TRUE)               
 }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
