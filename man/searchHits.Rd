\name{searchHits}
\alias{searchHits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Search for off targets in a sequence as DNAString
}
\description{
Search for off targets for given gRNAs, sequence and maximum mismatches
}
\usage{
searchHits(gRNAs, seqs, seqname, max.mismatch = 3, PAM.size = 3,
                 gRNA.size = 20, PAM = "NGG", PAM.pattern = "NNN$",
                 allowed.mismatch.PAM = 2, PAM.location = "3prime",
                 outfile)
}
\arguments{
  \item{gRNAs}{
DNAStringSet object containing a set of gRNAs. Please note the sequences 
must contain PAM appended after gRNAs, e.g., ATCGAAATTCGAGCCAATCCCGG where 
ATCGAAATTCGAGCCAATCC is the gRNA and CGG is the PAM
}
  \item{seqs}{
DNAString object containing a DNA sequence.
}
\item{seqname}{
Specify the name of the sequence 
}
  \item{max.mismatch}{
Maximum mismatch allowed in off target search, default 3. Warning: will be 
considerably slower if it is set to greater than 3 
}
  \item{PAM.size}{
Size of PAM, default 3
}
  \item{gRNA.size}{
Size of gRNA, default 20
}
\item{PAM}{
PAM as regular expression for appending to the gRNA, default NGG for SpCas9,
change to TTTN for cpf1.
}
\item{PAM.pattern}{
Regular expression of PAM, default N[A|G]G$ for spCas9.
For cpf1, ^TTTN since it is a 5 prime PAM sequence
}
 \item{allowed.mismatch.PAM}{
Maximum number of mismatches allowed in the offtargets comparing
to the PAM sequence. Default to 2 for NGG PAM
}
\item{PAM.location}{
PAM location relative to gRNA. For example, spCas9 PAM is located on the 3 prime while
cpf1 PAM is located on the 5 prime
}
\item{outfile}{
File path to temporarily store the search results
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
 if(interactive())
 {
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

    ### test PAM located at 5 prime
    all.gRNAs <- findgRNAs(inputFilePath = 
             system.file("extdata", "inputseq.fa", package = "CRISPRseek"),
             pairOutputFile = "pairedgRNAs.xls",
             findPairedgRNAOnly = FALSE,
             PAM = "TGT", PAM.location = "5prime")
     
    library("BSgenome.Hsapiens.UCSC.hg19")
         ### for speed reason, use max.mismatch = 0 for finding all targets with 
         ### all variants of PAM
    hits <- searchHits(all.gRNAs[1], BSgenomeName = Hsapiens, PAM.size = 3,
        max.mismatch = 0, chromToSearch = "chrX", PAM.location = "5prime",
        PAM = "^T[A|G]N", allowed.mismatch.PAM = 2)
    colnames(hits)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
