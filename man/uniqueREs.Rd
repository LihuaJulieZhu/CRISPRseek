\name{uniqueREs}
\alias{uniqueREs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Output restriction enzymes that recognize only the gRNA cleavage sites
}
\description{
For each identified gRNA, output restriction enzymes that recognize only the gRNA cleavage sites.
}
\usage{
uniqueREs(REcutDetails, summary, offTargets, scanUpstream = 100,
   scanDownstream = 100, BSgenomeName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{REcutDetails}{
REcutDetails stored in the REcutDetails.xls 
}
\item{summary}{
summary stored in the summary.xls}
\item{offTargets}{
offTargets stored in the offTargets.xls
}
\item{scanUpstream}{
upstream offset from the gRNA start, default 100
}
  \item{scanDownstream}{
downstream offset from the gRNA end, default 100
}

\item{BSgenomeName}{
BSgenome object. Please refer to available.genomes in BSgenome package. For 
example, BSgenome.Hsapiens.UCSC.hg19 for hg19, BSgenome.Mmusculus.UCSC.mm10 
for mm10, BSgenome.Celegans.UCSC.ce6 for ce6, BSgenome.Rnorvegicus.UCSC.rn5 
for rn5, BSgenome.Drerio.UCSC.danRer7 for Zv9, and 
BSgenome.Dmelanogaster.UCSC.dm3 for dm3
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 returns the RE sites that recognize only the gRNA cleavage sites for each gRNA.

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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
}
\examples{
    library("BSgenome.Hsapiens.UCSC.hg19")
    load(system.file("extdata", "ForTestinguniqueREs.RData",
            package = "CRISPRseek"))
    uniqueREs(results$REcutDetails, results$summary, results$offtarget,
	scanUpstream = 50,
        scanDownstream = 50, BSgenomeName = Hsapiens)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
