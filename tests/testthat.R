dependencies <- c("CRISPRseek", 
                  "BSgenome.Hsapiens.UCSC.hg19", 
                  "org.Hs.eg.db",
                  "TxDb.Hsapiens.UCSC.hg19.knownGene",
                  "BSgenome.Mmusculus.UCSC.mm10",
                  "org.Mm.eg.db",
                  "TxDb.Mmusculus.UCSC.mm10.knownGene",
                  "testthat")
sapply(dependencies, function(x) ifelse(do.call(require, list(x)), "okay", abort(x, "is not loaded for unit testing!")))
test_check("CRISPRseek")

