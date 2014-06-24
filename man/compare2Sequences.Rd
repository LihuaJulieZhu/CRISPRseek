\name{compare2Sequences}
\alias{compare2Sequences}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare 2 input sequences for possible guide RNAs (gRNAs)
}
\description{
Generate all possible guide RNAs (gRNAs) for two input sequences and generate 
    scores for potential off-targets in the other sequence.
}
\usage{
compare2Sequences(inputFile1Path, inputFile2Path, format = "fasta", 
    findgRNAsWithREcutOnly = FALSE, REpatternFile, minREpatternSize = 6,
    overlap.gRNA.positions = c(17, 18), findPairedgRNAOnly = FALSE, 
    min.gap = 0, max.gap = 20, gRNA.name.prefix = "gRNA", PAM.size = 3, 
    gRNA.size = 20, PAM = "NGG", PAM.pattern = "N[A|G]G$", max.mismatch =4,
    outputDir, 
    weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445,
    0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583), 
    overwrite = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputFile1Path}{
Sequence input file 1 path that contains one of the 
two sequences to be searched for potential gRNAs 
}
  \item{inputFile2Path}{
Sequence input file 2 path that contains one of the
two sequences to be searched for potential gRNAs 
}
  \item{format}{
Format of the input file, fasta and fastq are supported, default fasta
}
  \item{findgRNAsWithREcutOnly}{
Indicate whether to find gRNAs overlap with restriction enzyme recognition 
pattern
}
  \item{REpatternFile}{
File path containing restriction enzyme cut patters
}
  \item{minREpatternSize}{
Minimum restriction enzyme recognition pattern length required for the enzyme 
pattern to be searched for, default 6
}
  \item{overlap.gRNA.positions}{
The required overlap positions of gRNA and restriction enzyme cut site, 
default 17 and 18
}
  \item{findPairedgRNAOnly}{
Choose whether to only search for paired gRNAs in such an orientation that 
the first one is on minus strand called reverse gRNA and the second one is 
on plus strand called forward gRNA. TRUE or FALSE, default FALSE
}
  \item{min.gap}{
Minimum distance between two oppositely oriented gRNAs to be valid paired 
gRNAs. Default 0
}
  \item{max.gap}{
Maximum distance between two oppositely oriented gRNAs to be valid paired 
gRNAs. Default 20
}
  \item{gRNA.name.prefix}{
The prefix used when assign name to found gRNAs, default gRNA, short for 
guided RNA.
}
  \item{PAM.size}{
PAM length, default 3
}
  \item{gRNA.size}{
The size of the gRNA, default 20
}
  \item{PAM}{
PAM sequence after the gRNA, default NGG
}
  \item{PAM.pattern}{
Regular expression of PAM, default N[A|G]G$
}
  \item{max.mismatch}{
Maximum mismatch allowed to search the off targets in the other sequence,
default 4
}
  \item{outputDir}{
the directory where the sequence comparison results will be written to
}
  \item{weights}{
 numeric vector size of gRNA length, default c(0, 0, 0.014, 0, 0, 0.395, 
0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 
0.804, 0.685, 0.583) which is used in Hsu et al., 2013 cited in the reference 
section
}
  \item{overwrite}{
overwrite the existing files in the output directory or not, default TRUE 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Return a data frame with all potential gRNAs from both sequences.
In addition, a tab delimited file scoresFor2InputSequences.xls is also 
saved in the outputDir, sorted by scoreDiff descending.
\item{name}{name of the gRNA}
\item{gRNAPlusPAM}{gRNA plus PAM sequence}
\item{targetInSeq1}{target/off-target sequence in the 1st input sequence file}
\item{targetInSeq2}{target/off-target sequence in the 2nd input sequence file}
\item{scoreForSeq1}{score for the target sequence in the 1st input 
    sequence file}
\item{scoreForSeq2}{score for the target sequence in the 1st input 
    sequence file}
\item{mismatch.distance2PAM}{distances of mismatch to PAM, e.g.,
 14 means the mismatch is 14 bp away from PAM}
\item{n.mismatch}{number of mismatches between the off-target and the gRNA}
\item{targetSeqName}{the name of the input sequence where the target
    sequence is located}
\item{scoreDiff}{scoreForSeq1 - scoreForSeq2}
}
\references{
Patrick D Hsu, David A Scott, Joshua A Weinstein, F Ann Ran, Silvana 
Konermann, Vineeta Agarwala, Yinqing Li, Eli J Fine, Xuebing Wu, Ophir Shalem, 
Thomas J Cradick, Luciano A Marraffini, Gang Bao & Feng Zhang (2013) DNA 
targeting specificity of rNA-guided Cas9 nucleases. Nature Biotechnology 
31:827-834
}
\author{
Lihua Julie Zhu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
CRISPRseek
}
\examples{
    library(CRISPRseek)
    inputFile1Path <- system.file("extdata", "rs362331T.fa",
            package = "CRISPRseek")
    inputFile2Path <- system.file("extdata", "rs362331C.fa",
            package = "CRISPRseek")
    REpatternFile <- system.file("extdata", "NEBenzymes.fa", 
            package = "CRISPRseek")
    seqs <- compare2Sequences(inputFile1Path, inputFile2Path,
        outputDir = getwd(), 
        REpatternFile = REpatternFile, overwrite = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
