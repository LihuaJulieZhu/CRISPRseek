\name{offTargetAnalysisOfPeakRegions}
\alias{offTargetAnalysisOfPeakRegions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finding offtargets in given regions
}
\description{
Finding offtargets in given regions, such as peaks from GUIDE-seq
}
\usage{
offTargetAnalysisOfPeakRegions(gRNA, peaks, 
   format=c("fasta", "bed"),
   peaks.withHeader = FALSE, BSgenomeName,
   upstream = 50, downstream =50, PAM.size = 3, gRNA.size = 20,
   PAM = "NGG", PAM.pattern = "NNN$", max.mismatch = 6,
   outputDir, allowed.mismatch.PAM = 3, overwrite = TRUE,
   weights = c(0, 0, 0.014, 0, 0, 0.395,
0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615,
0.804, 0.685, 0.583)
   )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gRNA}{
gRNA input file path or a DNAStringSet object that contains gRNA plus PAM
sequences used for genome editing 
}
\item{peaks}{
peak input file path or a GenomicRanges object that contains genomic regions 
to be searched for potential offtargets
}
  \item{format}{
Format of the gRNA and peak input file. Currently, fasta and bed are supported
for gRNA and peak input file respectively
}
\item{peaks.withHeader}{
Indicate whether the peak input file contains header, default FALSE}
  \item{PAM.size}{
PAM length, default 3
}
  \item{gRNA.size}{
The size of the gRNA, default 20
}
  \item{PAM}{
PAM sequence after the gRNA, default NGG
}
  \item{BSgenomeName}{
BSgenome object. Please refer to available.genomes in BSgenome package. For 
example, BSgenome.Hsapiens.UCSC.hg19 for hg19, BSgenome.Mmusculus.UCSC.mm10 
for mm10, BSgenome.Celegans.UCSC.ce6 for ce6, BSgenome.Rnorvegicus.UCSC.rn5 
for rn5, BSgenome.Drerio.UCSC.danRer7 for Zv9, and 
BSgenome.Dmelanogaster.UCSC.dm3 for dm3
}
  \item{max.mismatch}{
Maximum mismatch allowed in off target search, default 6 
}
  \item{PAM.pattern}{
Regular expression of protospacer-adjacent motif (PAM), default N[A|G]G$
}
\item{allowed.mismatch.PAM}{
Number of degenerative bases in the PAM sequence, default to 2 for N[A|G]G PAM
}
  \item{outputDir}{
the directory where the off target analysis and reports will be written to
}
  \item{upstream}{
upstream offset from the peak start, default 50
}
  \item{downstream}{
downstream offset from the peak end, default 50
}
  \item{weights}{
a numeric vector size of gRNA length, default c(0, 0, 0.014, 0, 0, 0.395, 
0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 
0.804, 0.685, 0.583) which is used in Hsu et al., 2013 cited in the reference 
section
}
\item{overwrite}{
  overwrite the existing files in the output directory or not, default FALSE
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Same as compare2Sequences with one more tab-delimited file 
offTargetsInPeakRegions.xls containing all input peaks
with potential gRNA binding sites, mismatch number and positions, and 
predicted cleavage score.
}
\references{
Patrick D Hsu, David A Scott, Joshua A Weinstein, F Ann Ran, Silvana 
Konermann, Vineeta Agarwala, Yinqing Li, Eli J Fine, Xuebing Wu, Ophir Shalem, 
Thomas J Cradick, Luciano A Marraffini, Gang Bao & Feng Zhang (2013) DNA 
targeting specificity of rNA-guided Cas9 nucleases. Nature Biotechnology 
31:827-834
Lihua Julie Zhu, Benjamin R. Holmes, Neil Aronin and Michael Brodsky. CRISPRseek: a Bioconductor package to identify target-specific guide RNAs for CRISPR-Cas9 genome-editing systems. Plos One Sept 23rd 2014
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
   library("BSgenome.Hsapiens.UCSC.hg19")
   peaks <- system.file("extdata", "T2plus100OffTargets.bed",
      package = "CRISPRseek")
   gRNAs <- system.file("extdata", "T2.fa",
      package = "CRISPRseek")
   outputDir = getwd() 
   offTargets <- offTargetAnalysisOfPeakRegions(gRNA = gRNAs, peaks = peaks,
      format=c("fasta", "bed"),
      peaks.withHeader = TRUE, BSgenomeName = Hsapiens,
      upstream = 50, downstream =50, PAM.size = 3, gRNA.size = 20,
      PAM = "NGG", PAM.pattern = "NNN$", max.mismatch = 8,
      outputDir = outputDir,
      allowed.mismatch.PAM = 3, overwrite = TRUE
   )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
