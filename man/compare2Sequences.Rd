\name{compare2Sequences}
\alias{compare2Sequences}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare 2 input sequences/sequence sets for possible guide RNAs (gRNAs)
}
\description{
Generate all possible guide RNAs (gRNAs) for two input sequences, or two sets of sequences and generate scores for potential off-targets in the other sequence.
}
\usage{
compare2Sequences(inputFile1Path, inputFile2Path, 
    inputNames=c("Seq1", "Seq2"), 
    format = c("fasta", "fasta"), header=FALSE, findgRNAsWithREcutOnly = FALSE, 
    searchDirection=c("both","1to2", "2to1"), BSgenomeName,
    REpatternFile=system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek"),
    minREpatternSize = 6, findgRNAs = c(TRUE, TRUE), removegRNADetails = c(FALSE, FALSE),
    exportAllgRNAs = c("no", "all", "fasta", "genbank"), annotatePaired =  FALSE,
    overlap.gRNA.positions = c(17, 18), findPairedgRNAOnly = FALSE, 
    min.gap = 0, max.gap = 20, gRNA.name.prefix = "gRNA", PAM.size = 3, 
    gRNA.size = 20, PAM = "NGG", PAM.pattern = "N[A|G]G$",
    allowed.mismatch.PAM = 2, max.mismatch = 3,
    outputDir, upstream = 0, downstream = 0, 
    weights = c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445,
    0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583), 
    overwrite = FALSE, baseBeforegRNA = 4, 
    baseAfterPAM = 3, featureWeightMatrixFile = system.file("extdata", 
       "DoenchNBT2014.csv", package = "CRISPRseek"), foldgRNAs = FALSE, 
    gRNA.backbone="GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU",
    temperature = 37,
    scoring.method = c("Hsu-Zhang", "CFDscore"),
        subPAM.activity = hash( AA =0,
          AC =   0,
          AG = 0.259259259,
          AT = 0,
          CA = 0,
          CC = 0,
          CG = 0.107142857,
          CT = 0,
          GA = 0.069444444,
          GC = 0.022222222,
          GG = 1,
          GT = 0.016129032,
          TA = 0,
          TC = 0,
          TG = 0.038961039,
          TT = 0),
     subPAM.position = c(22, 23),
     mismatch.activity.file = system.file("extdata",
         "NatureBiot2016SuppTable19DoenchRoot.csv",
         package = "CRISPRseek")    
    )
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
  \item{inputNames}{
Name of the input sequences when inputFile1Path and inputFile2Path are 
DNAStringSet instead of file path
}
  \item{format}{
Format of the input files, fasta, fastq and bed format are supported, default fasta
}
  \item{header}{
Indicate whether the input file contains header, default FALSE, only applies to bed format}
  \item{findgRNAsWithREcutOnly}{
Indicate whether to find gRNAs overlap with restriction enzyme recognition 
pattern
}
  \item{searchDirection}{
Indicate whether perfrom gRNA in both sequences and off-target search against each other (both)
or search gRNA in input1 and off-target analysis in input2 (1to2), or vice versa (2to1)
}
\item{BSgenomeName}{
BSgenome object. Please refer to available.genomes in BSgenome package. For 
example, BSgenome.Hsapiens.UCSC.hg19 for hg19, BSgenome.Mmusculus.UCSC.mm10 
for mm10, BSgenome.Celegans.UCSC.ce6 for ce6, BSgenome.Rnorvegicus.UCSC.rn5 
for rn5, BSgenome.Drerio.UCSC.danRer7 for Zv9, and 
BSgenome.Dmelanogaster.UCSC.dm3 for dm3
}
  \item{REpatternFile}{
File path containing restriction enzyme cut patters
}
  \item{minREpatternSize}{
Minimum restriction enzyme recognition pattern length required for the enzyme 
pattern to be searched for, default 6
}
  \item{findgRNAs}{
Indicate whether to find gRNAs from the sequences in the input file or skip 
the step of finding gRNAs, default TRUE for both input sequences. 
Set it to FALSE if the input file contains user selected gRNAs plus PAM already. 
}
\item{removegRNADetails}{
Indicate whether to remove the detailed gRNA information such as efficacy file and 
restriction enzyme cut sites, default false for both input sequences.
Set it to TRUE if the input file contains the user selected gRNAs plus PAM already.
}
\item{exportAllgRNAs}{
Indicate whether to output all potential gRNAs to a file in fasta format,
genbank format or both. Default to no.
}
\item{annotatePaired}{Indicate whether to output paired information, default to FALSE}
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
\item{allowed.mismatch.PAM}{
Number of degenerative bases in the PAM sequence, default to 2 for N[A|G]G PAM
}
  \item{max.mismatch}{
Maximum mismatch allowed to search the off targets in the other sequence,
default 3 
}
  \item{outputDir}{
the directory where the sequence comparison results will be written to
}
  \item{upstream}{
upstream offset from the bed input starts to search for gRNA and/or offtargets, default 0
}
  \item{downstream}{
downstream offset from the bed input ends to search for gRNA and/or offtargets, default 0 
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
\item{baseBeforegRNA}{
Number of bases before gRNA used for calculating gRNA efficiency, default 4
}
\item{baseAfterPAM}{
Number of bases after PAM used for calculating gRNA efficiency, default 3
}
\item{featureWeightMatrixFile}{
Feature weight matrix file used for calculating gRNA efficiency. By default DoenchNBT2014 weight matrix is used. To use alternative weight matrix file, please input a csv file with first column containing significant features and the second column containing the corresponding weights for the features. Please see Doench et al., 2014 for details.
}
\item{foldgRNAs}{
Default FALSE. If set to TRUE, summary file will contain minimum free energy of the secondary
structure of gRNA with gRNA backbone from GeneRfold package provided that GeneRfold package has been installed.
}
\item{gRNA.backbone}{
gRNA backbone constant region sequence. Default to the
sequence in Sp gRNA backbone.
}
  \item{temperature}{
temperature in celsius. Default to 37 celsius. 
}
\item{scoring.method}{
   Indicates which method to use for offtarget cleavage rate estimation, currently
two methods are supported,  Hsu-Zhang and CFDscore
}
\item{subPAM.activity}{
   Applicable only when scoring.method is set to CFDscore
   A hash to represent the cleavage rate for each alternative sub PAM sequence relative
to preferred PAM sequence
}
\item{subPAM.position}{
   Applicable only when scoring.method is set to CFDscore
The start and end positions of the sub PAM. Default to 22 and 23 for SP with 20bp gRNA and NGG as preferred PAM
}
\item{mismatch.activity.file}{
Applicable only when scoring.method is set to CFDscore
A comma separated (csv) file containing the cleavage rates for all possible types of
single nucleotide mismatche at each position of the gRNA. By default, using the supplemental Table 19 from Doench et al., Nature Biotechnology 2016
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
\item{targetInSeq1}{target/off-target sequence including PAM in the 1st 
	input sequence file}
\item{targetInSeq2}{target/off-target sequence incuding PAM in the 2nd 
	input sequence file}
\item{guideAlignment2Offtarget}{alignment of gRNA to the other input 
sequence (off-target sequence)}
\item{offTargetStrand}{strand of the other sequence (off-target sequence)
	the gRNA align to}
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
\item{bracket.notation}{folded gRNA in bracket notation}
\item{mfe.sgRNA}{minimum free energy of sgRNA}
\item{mfe.diff}{mfe.sgRNA-mfe.backbone}
\item{mfe.backbone}{minimum free energy of the gRNA backbone by itself}
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
