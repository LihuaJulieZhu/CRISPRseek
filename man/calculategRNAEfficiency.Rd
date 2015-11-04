\name{calculategRNAEfficiency}
\alias{calculategRNAEfficiency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Calculate gRNA Efficiency
}
\description{
    Calculate gRNA Efficiency for a given set of sequences and feature weight matrix
}
\usage{
calculategRNAEfficiency(extendedSequence, 
    baseBeforegRNA, featureWeightMatrix, gRNA.size = 20,
    enable.multicore = FALSE, n.cores.max = 6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{extendedSequence}{
 Sequences containing gRNA plus PAM plus flanking sequences. Each sequence should be long enough for building features specified in the featureWeightMatrix
}
\item{baseBeforegRNA}{
Number of bases before gRNA used for calculating gRNA efficiency, default 4
}
  \item{featureWeightMatrix}{
  a data frame with the first column containing significant features and the second column containing the weight of corresponding features. In the following example, DoenchNBT2014 weight matrix is used. Briefly, features include INTERCEPT,GC_LOW (penalty for low GC content in the gRNA sequence), GC_HIGH (penalty for high GC content in the gRNA sequence), G02 (means G at second position of the extendedSequence), GT02 (means GT di-nucleotides starting at 2nd position of the extendedSequence). To understand how is the feature weight matrix is identified, or how to use alternative feature weight matrix file, please see Doench et al., 2014 for details. 
}
    \item{gRNA.size}{
The size of the gRNA, default 20
}
\item{enable.multicore}{Indicate whether enable parallel processing, default FALSE.
For super long sequences with lots of gRNAs, suggest set it to TRUE}
\item{n.cores.max}{
Indicating maximum number of cores to use in multi core mode, 
    i.e., parallel processing, default 6. Please set it to 1 to disable
    multicore processing for small dataset.
}
 }

\value{
DNAStringSet consists of potential gRNAs that can be input to filtergRNAs 
function directly
}
\references{
Doench JG, Hartenian E, Graham DB, Tothova Z, Hegde M, Smith I, Sullender M, Ebert BL, Xavier RJ, Root DE. Rational design of highly active sgRNAs for CRISPR-Cas9-mediated gene inactivation. Nat Biotechnol. 2014 Sep 3. doi: 10.1038 nbt.3026
http://www.broadinstitute.org/rnai/public/analysis-tools/sgrna-design
}
\author{
Lihua Julie Zhu
}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
offTargetAnalysis
}
\examples{
 	extendedSequence <- c("TGGATTGTATAATCAGCATGGATTTGGAAC",
		"TCAACGAGGATATTCTCAGGCTTCAGGTCC",
		"GTTACCTGAATTTGACCTGCTCGGAGGTAA",
		"CTTGGTGTGGCTTCCTTTAAGACATGGAGC",
		"CATACAGGCATTGAAGAAGAATTTAGGCCT",
		"AGTACTATACATTTGGCTTAGATTTGGCGG",
		"TTTTCCAGATAGCCGATCTTGGTGTGGCTT",
		"AAGAAGGGAACTATTCGCTGGTGATGGAGT"
	)
	featureWeightMatrixFile <- system.file("extdata", "DoenchNBT2014.csv", 
		package = "CRISPRseek")
	featureWeightMatrix <- read.csv(featureWeightMatrixFile, header=TRUE)
	calculategRNAEfficiency(extendedSequence, baseBeforegRNA = 4, 
		featureWeightMatrix, gRNA.size = 20)

 }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
