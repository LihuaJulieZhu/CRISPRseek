\name{writeHits2}
\alias{writeHits2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write the hits of sequence search to a file
}
\description{
write the hits of sequence search to a file, internal function used by 
searchHits
}
\usage{
writeHits2(gRNA, seqname, matches, strand, file, gRNA.size = 20, 
    PAM = "NGG", PAM.pattern = "N[A|G]G$",
    max.mismatch = 4, chrom.len, append = FALSE,
    PAM.location = "3prime", PAM.size = 3,
    allowed.mismatch.PAM = 1L,
    BSgenomeName, baseEditing = FALSE, targetBase = "C", editingWindow = 5:13)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gRNA}{
DNAString object with gRNA sequence with PAM appended immediately after,e.g.,
ACGTACGTACGTACTGACGTCGG with 20bp gRNA sequence plus 3bp PAM sequence CGG
}
  \item{seqname}{
chromosome name as character, e.g., chr1
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
PAM as regular expression for filtering the hits, default NGG for spCas9.
For cpf1, TTTN.
}
  \item{PAM.pattern}{
Regular expression of protospacer-adjacent motif (PAM), default N[A|G]G$ for spCas9.
For cpf1, ^TTTN since it is a 5 prime PAM sequence
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
Number of degenerative bases in the PAM sequence, default to 1 for N[A|G]G PAM
}
 \item{BSgenomeName}{
BSgenome object. Please refer to available.genomes in BSgenome package. For
example, BSgenome.Hsapiens.UCSC.hg19 for hg19, BSgenome.Mmusculus.UCSC.mm10
for mm10, BSgenome.Celegans.UCSC.ce6 for ce6, BSgenome.Rnorvegicus.UCSC.rn5
for rn5, BSgenome.Drerio.UCSC.danRer7 for Zv9, and
BSgenome.Dmelanogaster.UCSC.dm3 for dm3
}
\item{baseEditing}{
Indicate whether to design gRNAs for base editing. Default to FALSE
If TRUE, please set baseEditing = TRUE, targetBase and editingWidow accordingly.
}
\item{targetBase}{
Applicable only when baseEditing is set to TRUE. It is used to indicate the target
base for base editing systems, default to C for converting C to T in the CBE system.
Please change it to A if you intend to use the ABE system.
}
\item{editingWindow}{
Applicable only when baseEditing is set to TRUE. It is used to indicate the
effective editing window to consider for the offtargets search only, 
default to 5 to 13 which is for the original CBE system.
Please change it accordingly if the system you use have a different editing window,
or you would like to include offtargets with the target base in a larger editing window.
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
    library("BSgenome.Hsapiens.UCSC.hg19")
    gRNAPlusPAM <- DNAString("ACGTACGTACGTACTGACGTCGG")
    x <- DNAString("AAGCGCGATATGACGTACGTACGTACTGACGTCGG")
    chrom.len <- nchar(as.character(x))
    m <- matchPattern(gRNAPlusPAM, x)
    names(m) <- "testing"
    writeHits2(gRNA = gRNAPlusPAM, seqname = "chr1", 
        PAM = "NGG", PAM.pattern = "NNN$", allowed.mismatch.PAM = 2,
        matches = m, strand = "+", file = "exampleWriteHits.txt", 
        chrom.len = chrom.len, append = FALSE, BSgenomeName = Hsapiens)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
