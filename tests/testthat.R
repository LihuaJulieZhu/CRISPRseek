Sys.setenv(R_TESTS="")
require("CRISPRseek") || stop("unable to load Package:CRISPRseek")
require("BSgenome.Hsapiens.UCSC.hg19") || 
  stop("unable to load Package:BSgenome.Hsapiens.UCSC.hg19")
require("org.Hs.eg.db") || stop("unable to load Package:org.Hs.eg.db")
require("TxDb.Hsapiens.UCSC.hg19.knownGene") || 
  stop("unable to load TxDb.Hsapiens.UCSC.hg19.knownGene")
require("BSgenome.Mmusculus.UCSC.mm10") || 
  stop("unable to load Package: BSgenome.Mmusculus.UCSC.mm10")
require("TxDb.Mmusculus.UCSC.mm10.knownGene") || 
  stop("unable to load Package:TxDb.Mmusculus.UCSC.mm10.knownGene")
require("org.Mm.eg.db") || 
  stop("unable to load Package:org.Mm.eg.db")
require("testthat") || stop("unable to load testthat")
test_check("CRISPRseek")
