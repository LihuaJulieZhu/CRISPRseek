\name{offTargetAnalysis}
\alias{offTargetAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Design of target-specific guide RNAs for CRISPR-Cas9 system in one function
}
\description{
Design of target-specific guide RNAs (gRNAs) for CRISPR-Cas9 system by automatically 
calling findgRNAs, filtergRNAs, searchHits, buildFeatureVectorForScoring, 
getOfftargetScore, filterOfftarget and generate reports.
}
\usage{
offTargetAnalysis(inputFilePath, format = "fasta", findgRNAs = TRUE,
    exportAllgRNAs = c("all", "fasta", "genbank", "no"),
    findgRNAsWithREcutOnly = TRUE, REpatternFile, minREpatternSize = 6, 
    overlap.gRNA.positions = c(17, 18), findPairedgRNAOnly = TRUE, 
    min.gap = 0, max.gap = 20, gRNA.name.prefix = "gRNA", PAM.size = 3, 
    gRNA.size = 20,  PAM = "NGG",  BSgenomeName, chromToSearch = "all", 
    max.mismatch = 4, PAM.pattern = "N[A|G]G$", gRNA.pattern = "",
    min.score = 0.5, topN = 100, 
    topN.OfftargetTotalScore = 10, annotateExon = TRUE,
    txdb, outputDir, fetchSequence = TRUE, upstream = 200, downstream = 200, 
    weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 
    0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583),
    overwrite = FALSE)
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
  \item{findgRNAs}{
Indicate whether to find gRNAs from the sequences in the input file or skip 
the step of finding gRNAs, default TRUE. Set it to FALSE if the input file 
contains user selected gRNAs plus PAM already. 
}
\item{exportAllgRNAs}{
Indicate whether to output all potential gRNAs to a file in fasta format,
genbank format or both. Default to both.
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
  \item{BSgenomeName}{
BSgenome object. Please refer to available.genomes in BSgenome package. For 
example, BSgenome.Hsapiens.UCSC.hg19 for hg19, BSgenome.Mmusculus.UCSC.mm10 
for mm10, BSgenome.Celegans.UCSC.ce6 for ce6, BSgenome.Rnorvegicus.UCSC.rn5 
for rn5, BSgenome.Drerio.UCSC.danRer7 for Zv9, and 
BSgenome.Dmelanogaster.UCSC.dm3 for dm3
}
\item{chromToSearch}{
Specify the chromosome to search, default to all, meaning search all 
chromosomes. For example, chrX indicates searching for matching in chromosome 
X only
}
  \item{max.mismatch}{
Maximum mismatch allowed in off target search, default 4. Warning: will be 
considerably slower if set >4
}
  \item{PAM.pattern}{
Regular expression of protospacer-adjacent motif (PAM), default N[A|G]G$
}
\item{gRNA.pattern}{
Regular expression or IUPAC Extended Genetic Alphabet to represent gRNA 
pattern, default is no restriction. To specify that the gRNA must start 
with GG for example, then set it to ^GG. Please see help(translatePattern) for
a list of IUPAC Extended Genetic Alphabet.
}
  \item{min.score}{
minimum score of an off target to included in the final output, default 0.5
}
  \item{topN}{
top N off targets to be included in the final output, default 100
}
  \item{topN.OfftargetTotalScore}{
top N off target used to calculate the total off target score, default 10
}
  \item{annotateExon}{
Choose whether or not to indicate whether the off target is inside an exon or 
not, default TRUE
}
  \item{txdb}{
TranscriptDb object, for creating and using TranscriptDb object, please refer 
to GenomicFeatures package. For a list of existing TranscriptDb object, please 
search for annotation package starting with Txdb at 
http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData, 
such as TxDb.Rnorvegicus.UCSC.rn5.refGene for rat, 
TxDb.Mmusculus.UCSC.mm10.knownGene for mouse, 
TxDb.Hsapiens.UCSC.hg19.knownGene for human, 
TxDb.Dmelanogaster.UCSC.dm3.ensGene for Drosophila and 
TxDb.Celegans.UCSC.ce6.ensGene for C.elegans

}
  \item{outputDir}{
the directory where the off target analysis and reports will be written to
}
  \item{fetchSequence}{
Fetch flank sequence of off target or not, default TRUE
}
  \item{upstream}{
upstream offset from the off target start, default 200
}
  \item{downstream}{
downstream offset from the off target end, default 200
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
Four tab delimited files are generated in the output directory: 
OfftargetAnalysis.xls (detailed information of off targets), 
Summary.xls (summary of the gRNAs), REcutDetails.xls 
(restriction enzyme cut sites of each gRNA), and pairedgRNAs.xls 
(potential paired gRNAs)
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
	library("BSgenome.Hsapiens.UCSC.hg19")
	library(TxDb.Hsapiens.UCSC.hg19.knownGene)
	outputDir <- getwd()
	inputFilePath <- system.file("extdata", "inputseq.fa",
            package = "CRISPRseek")
	REpatternFile <- system.file("extdata", "NEBenzymes.fa", 
            package = "CRISPRseek")
	offTargetAnalysis(inputFilePath, findgRNAsWithREcutOnly = TRUE, 
            REpatternFile = REpatternFile, findPairedgRNAOnly = FALSE, 
            BSgenomeName = Hsapiens, chromToSearch = "chrX",
            txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, max.mismatch = 1, 
            outputDir = outputDir, overwrite = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
