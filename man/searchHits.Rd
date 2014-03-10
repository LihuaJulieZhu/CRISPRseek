\name{searchHits}
\alias{searchHits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Search for off targets
}
\description{
Search for off targets for given gRNAs, BSgenome and maximum mismatches
}
\usage{
searchHits(gRNAs, BSgenomeName, chromToSearch = "all", max.mismatch = 4, 
    PAM.size = 3, gRNA.size = 20, PAM = "N[A|G]G$")
}
\arguments{
  \item{gRNAs}{
DNAStringSet object containing a set of gRNAs. Please note the sequences 
must contain PAM appended after gRNAs, e.g., ATCGAAATTCGAGCCAATCCCGG where 
ATCGAAATTCGAGCCAATCC is the gRNA and CGG is the PAM
}
  \item{BSgenomeName}{
BSgenome object. Please refer to available.genomes in BSgenome package. 
For example, BSgenome.Hsapiens.UCSC.hg19 for hg19, BSgenome.Mmusculus.UCSC.mm10 
for mm10, BSgenome.Celegans.UCSC.ce6 for ce6, BSgenome.Rnorvegicus.UCSC.rn5 
for rn5, and BSgenome.Dmelanogaster.UCSC.dm3 for dm3
}
	\item{chromToSearch}{
Specify the chromosome to search, default to all, meaning search all 
chromosomes. For example, chrX indicates searching for matching in 
chromosome X only
}
  \item{max.mismatch}{
Maximum mismatch allowed in off target search, default 4. Warning: will be 
considerably slower if it is set to greater than 4
}
  \item{PAM.size}{
Size of PAM, default 3
}
  \item{gRNA.size}{
Size of gRNA, default 20
}
  \item{PAM}{
Regular expression of PAM, default N[A|G]G$
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a data frame contains IsMismatch.posX (indicator variable indicating whether 
this position X is mismatch or not, 1 means yes and 0 means not, 
X = 1 to gRNA.size) representing all positions in the gRNA),strand 
(strand of the match, + for plus and - for minus strand),	chrom 
(chromosome of the off target), chromStart (start position of the off target),	
chromEnd (end position of the off target),name (gRNA name), 
gRNAPlusPAM (gRNA sequence with PAM sequence concatenated), 
OffTargetSequence (the genomic sequence of the off target), n.mismatch 
(number of mismatches between the off target and the gRNA), forViewInUCSC 
(string for viewing in UCSC genome browser, e.g., chr14:31665685-31665707),
score (set to 100, and will be updated in getOfftargetScore)
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
\seealso{
offTargetAnalysis
}
\examples{
    all.gRNAs <- findgRNAs(inputFilePath = 
        system.file("extdata", "inputseq.fa", package = "CRISPRseek"),
        pairOutputFile = "pairedgRNAs.xls",
	findPairedgRNAOnly = TRUE)

    library("BSgenome.Hsapiens.UCSC.hg19")
    ### for speed reason, use max.mismatch = 0 for finding all targets with 
    ### all variants of PAM
    hits <- searchHits(all.gRNAs[1], BSgenomeName = Hsapiens,
        max.mismatch = 0, chromToSearch = "chrX")
    colnames(hits)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
