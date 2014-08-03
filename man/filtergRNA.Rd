\name{filtergRNAs}
\alias{filtergRNAs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter gRNAs
}
\description{
Filter gRNAs containing restriction enzyme cut site
}
\usage{
filtergRNAs(all.gRNAs, pairOutputFile = "", 
    findgRNAsWithREcutOnly = FALSE,
    REpatternFile = system.file("extdata", "NEBenzymes.fa",
        package = "CRISPRseek"), format = "fasta", 
    minREpatternSize = 6, overlap.gRNA.positions = c(17, 18))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{all.gRNAs}{
gRNAs as DNAStringSet, such as the output from findgRNAs
}
  \item{pairOutputFile}{
File path with paired gRNAs
}
  \item{findgRNAsWithREcutOnly}{
Indicate whether to find gRNAs overlap with restriction enzyme recognition 
pattern
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
  \item{overlap.gRNA.positions}{
The required overlap positions of gRNA and restriction enzyme cut site, 
default 17 and 18
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{gRNAs.withRE }{gRNAs as DNAStringSet that passed the filter criteria}
 \item{gRNAREcutDetails}{a data frame that contains a set of gRNAs 
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
    all.gRNAs <- findgRNAs(
        inputFilePath = system.file("extdata", "inputseq.fa", 
        package = "CRISPRseek"),
        pairOutputFile = "testpairedgRNAs.xls",
        findPairedgRNAOnly = TRUE)

    gRNAs.RE <- filtergRNAs(all.gRNAs = all.gRNAs,
        pairOutputFile = "testpairedgRNAs.xls",
        REpatternFile = system.file("extdata", "NEBenzymes.fa", 
        package = "CRISPRseek"))

    gRNAs  <- gRNAs.RE$gRNAs.withRE
    restriction.enzyme.cut.sites <- gRNAs.RE$gRNAREcutDetails
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
