\name{writeHits}
\alias{writeHits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write the hits of sequence search from a sequence to a file
}
\description{
write the hits of sequence search from a sequence instead of BSgenome to a file,
internal function used by 
searchHits
}
\usage{
writeHits(gRNA, seqname, matches, strand, file, gRNA.size = 20L,
     PAM = "NGG", PAM.pattern = "N[A|G]G$", max.mismatch = 4L,
     chrom.len, append = FALSE, PAM.location = "3prime", 
     PAM.size = 3L, allowed.mismatch.PAM = 1L,
     seqs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gRNA}{
DNAString object with gRNA sequence with PAM appended immediately after,e.g.,
ACGTACGTACGTACTGACGTCGG with 20bp gRNA sequence plus 3bp PAM sequence CGG
}
  \item{seqname}{
sequence name as character 
}
  \item{matches}{
XStringViews object storing matched chromosome locations
}
  \item{strand}{
strand of the match, + for plus strand and - for minus strand
}
  \item{file}{
file path where the hits is written to
}
  \item{gRNA.size}{
gRNA size, default 20
}
\item{PAM}{
PAM as regular expression for appending to the gRNA, default NGG for SpCas9,
change to TTTN for cpf1.
}
  \item{PAM.pattern}{
PAM as regular expression for filtering the hits, default N[A|G]G$ for spCas9.
For cpf1, ^TTTN since it is a 5 prime PAM sequence.
}
  \item{max.mismatch}{
maximum mismatch allowed within the gRNA (excluding PAM portion) for 
filtering the hits, default 4
}
  \item{chrom.len}{
length of the matched chromosome
}
  \item{append}{
TRUE if append to existing file, false if start a new file
}
\item{PAM.location}{
PAM location relative to gRNA. For example, spCas9 PAM is located on the 3 prime while
cpf1 PAM is located on the 5 prime
}
\item{PAM.size}{
Size of PAM, default 3
}
\item{allowed.mismatch.PAM}{
Maximum number of mismatches allowed in the offtargets comparing
to the PAM sequence. Default to 1 for NGG PAM
}
\item{seqs}{
DNAString object containing a DNA sequence.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
results are saved in the file specified by file
}
\references{
http://bioconductor.org/packages/2.8/bioc/vignettes/BSgenome/inst/doc/
GenomeSearching.pdf}
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
 if(interactive())
 {
    gRNAPlusPAM <- DNAString("ACGTACGTACGTACTGACGTCGG")
    x <- DNAString("AAGCGCGATATGACGTACGTACGTACTGACGTCGG")
    chrom.len <- nchar(as.character(x))
    m <- matchPattern(gRNAPlusPAM, x)
    names(m) <- "testing"
    writeHits(gRNA = gRNAPlusPAM, seqname = "chr1", 
        matches = m, strand = "+", file = "exampleWriteHits.txt", 
        chrom.len = chrom.len, append = FALSE)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
