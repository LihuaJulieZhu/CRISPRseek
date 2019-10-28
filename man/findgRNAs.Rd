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
findgRNAs(inputFilePath,  baseEditing = FALSE, targetBase = "C", editingWindow = 5:13,
    format = "fasta", PAM = "NGG", PAM.size = 3, 
    findPairedgRNAOnly = FALSE, annotatePaired = TRUE,
    paired.orientation = c("PAMout", "PAMin"), enable.multicore = FALSE,
    n.cores.max = 6, gRNA.pattern = "", gRNA.size = 20, 
    overlap.gRNA.positions = c(17,18), min.gap = 0, max.gap = 20,
    pairOutputFile, name.prefix = "", featureWeightMatrixFile = system.file("extdata", 
        "DoenchNBT2014.csv", package = "CRISPRseek"), baseBeforegRNA = 4, 
        baseAfterPAM = 3, calculategRNAEfficacy = FALSE, efficacyFile,
    PAM.location = "3prime", rule.set = c("Root_RuleSet1_2014", "Root_RuleSet2_2016"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputFilePath}{
 Sequence input file path or a DNAStringSet object that contains sequences 
 to be searched for potential gRNAs
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
effective editing window, default to 5 to 13 which is for the original CBE system.
Please change it accordingly if the system you use have a different editing window.
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
\item{annotatePaired}{Indicate whether to output paired information, default TRUE}
\item{paired.orientation}{PAMin orientation means the two adjacent PAMs on the sense and 
antisense strands face inwards towards each other like N21GG and CCN21 whereas PAMout orientation
means they face away from each other like CCN21 and N21GG}
\item{enable.multicore}{Indicate whether enable parallel processing, default FALSE.
For super long sequences with lots of gRNAs, suggest set it to TRUE}
\item{n.cores.max}{
Indicating maximum number of cores to use in multi core mode, 
    i.e., parallel processing, default 6. Please set it to 1 to disable
    multicore processing for small dataset.
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
\item{overlap.gRNA.positions}{
The required overlap positions of gRNA and restriction enzyme cut site, 
default 17 and 18
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
 \item{baseBeforegRNA}{
Number of bases before gRNA used for calculating gRNA efficiency, default 4 for spCas9
Please note, for PAM located on the 5 prime, need to specify the number of bases before the PAM
sequence plus PAM size.
}
\item{baseAfterPAM}{
Number of bases after PAM used for calculating gRNA efficiency, default 3 for spCas9
Please note, for PAM located on the 5 prime, need to include the length of the gRNA plus
the extended sequence on the 3 prime
}
\item{featureWeightMatrixFile}{
Feature weight matrix file used for calculating gRNA efficiency. By default DoenchNBT2014 weight matrix is used. To use alternative weight matrix file, please input a csv file with first column containing significant features and the second column containing the corresponding weights for the features. Please see Doench et al., 2014 for details.
}
\item{calculategRNAEfficacy}{
Default to FALSE, not to calculate gRNA efficacy
}
\item{efficacyFile}{
File path to write gRNA efficacies
}
\item{PAM.location}{
PAM location relative to gRNA. For example, spCas9 PAM is located on the 3 prime while
cpf1 PAM is located on the 5 prime
}
\item{rule.set}{
Specify a rule set scoring system for calculating gRNA efficacy.
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
    findgRNAs(inputFilePath = system.file("extdata", 
        "cpf1.fa", package = "CRISPRseek"), 
        findPairedgRNAOnly=FALSE, 
        pairOutputFile = "testpairedgRNAs-cpf1.xls", 
        PAM="TTTN", PAM.location = "5prime", PAM.size = 4, 
        overlap.gRNA.positions = c(19,23),
        baseBeforegRNA = 8, baseAfterPAM = 23,
        calculategRNAEfficacy= TRUE, 
        featureWeightMatrixFile = system.file("extdata", 
            "DoenchNBT2014.csv", package = "CRISPRseek"),
       efficacyFile = "testcpf1Efficacy.xls")

    findgRNAs(inputFilePath = system.file("extdata", 
             "cpf1.fa", package = "CRISPRseek"), 
             findPairedgRNAOnly=FALSE, 
             pairOutputFile = "testpairedgRNAs-cpf1.xls", 
             PAM="TTTN", PAM.location = "5prime", PAM.size = 4, 
             overlap.gRNA.positions = c(19,23),
             baseBeforegRNA = 8, baseAfterPAM = 23,
             calculategRNAEfficacy= TRUE, 
             featureWeightMatrixFile = system.file("extdata", 
                 "DoenchNBT2014.csv", package = "CRISPRseek"),
            efficacyFile = "testcpf1Efficacy.xls", baseEditing =  TRUE, 
            editingWindow=20, targetBase = "X")

    findgRNAs(inputFilePath = system.file("extdata", 
             "cpf1.fa", package = "CRISPRseek"), 
             findPairedgRNAOnly=FALSE, 
             pairOutputFile = "testpairedgRNAs-cpf1.xls", 
             PAM="TTTN", PAM.location = "5prime", PAM.size = 4, 
             overlap.gRNA.positions = c(19,23),
             baseBeforegRNA = 8, baseAfterPAM = 23,
             calculategRNAEfficacy= TRUE, 
             featureWeightMatrixFile = system.file("extdata", 
                 "DoenchNBT2014.csv", package = "CRISPRseek"),
            efficacyFile = "testcpf1Efficacy.xls", baseEditing =  TRUE,
            editingWindow=20, targetBase = "C")
 }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
