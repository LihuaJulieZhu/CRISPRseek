\name{findgRNAs}
\alias{findgRNAs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Find potential gRNAs
}
\description{
	Find potential gRNAs for an input file containing sequences in fasta 
    format
}
\usage{
findgRNAs(inputFilePath, format = "fasta", PAM = "NGG", PAM.size = 3, 
    findPairedgRNAOnly = FALSE, gRNA.pattern = "", gRNA.size = 20, min.gap = 0, max.gap = 20,
    pairOutputFile, name.prefix = "gRNA")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputFilePath}{
 Sequence input file path or a DNAStringSet object that contains sequences 
 to be searched for potential gRNAs
}
  \item{format}{
Format of the input file, fasta and fastq are supported, default fasta
}
  \item{PAM}{
protospacer-adjacent motif (PAM) sequence after the gRNA, default NGG
}
  \item{PAM.size}{
PAM length, default 3
}
  \item{findPairedgRNAOnly}{
Choose whether to only search for paired gRNAs in such an orientation that 
the first one is on minus strand called reverse gRNA and the second one is 
on plus strand called forward gRNA. TRUE or FALSE, default FALSE
}
\item{gRNA.pattern}{
Regular expression or IUPAC Extended Genetic Alphabet to represent gRNA 
pattern, default is no restriction. To specify that the gRNA must start 
with GG for example, then set it to ^GG. Please see help(translatePattern) for
a list of IUPAC Extended Genetic Alphabet.
}
  \item{gRNA.size}{
The size of the gRNA, default 20
}
  \item{min.gap}{
Minimum distance between two oppositely oriented gRNAs to be valid paired 
gRNAs. Default 0
}
  \item{max.gap}{
Maximum distance between two oppositely oriented gRNAs to be valid paired 
gRNAs. Default 20
}
  \item{pairOutputFile}{
The output file for writing paired gRNA information to
}
  \item{name.prefix}{
The prefix used when assign name to found gRNAs, default gRNA, short for
guided RNA.
}
}
\details{
If users already has a fasta file that contains a set of potential gRNAs, 
then users can call filergRNAs directly although the easiest way is to call 
the one-stop-shopping function OffTargetAnalysis with findgRNAs set to FALSE.
}
\value{
DNAStringSet consists of potential gRNAs that can be input to filtergRNAs 
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
    findgRNAs(inputFilePath = system.file("extdata",
        "inputseq.fa", package = "CRISPRseek"),
        pairOutputFile = "testpairedgRNAs.xls",
        findPairedgRNAOnly = TRUE)               
 }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
