test_need_tensorflow <- function() {
  test_that("test_findgRNAsCpf1", {
    pairOutputFile <- tempfile()
    efficacyFile <- tempfile()
    gRNAs <- findgRNAs(inputFilePath = system.file("extdata", "cpf1.fa", package = "CRISPRseek"),
                       findPairedgRNAOnly = FALSE,
                       pairOutputFile = pairOutputFile, 
                       PAM = "TTTN",
                       PAM.location = "5prime", 
                       PAM.size = 4, 
                       overlap.gRNA.positions = c(19, 23), 
                       baseBeforegRNA = 8, 
                       baseAfterPAM = 26, 
                       calculategRNAEfficacy = TRUE,
                       rule.set = "DeepCpf1", 
                       efficacyFile = efficacyFile)
    gRNAs_truth <- readRDS(system.file("extdata", "gRNAsCpf1.RDS", package = "CRISPRseek"))
    expect_equal(as.character(gRNAs), as.character(gRNAs_truth))
    
    effi_truth <- readRDS(system.file("extdata", "Cpf1Eff.RDS", package = "CRISPRseek"))
    effi <- read.table(efficacyFile, sep = "\t", header = TRUE)
    expect_equal(effi_truth, effi, tolerance = 0.001)
  })
}

# test_need_tensorflow()

test_that("test_findgRNAsSpCas9", {
  pairOutputFile <- tempfile()
  efficacyFile <- file.path(tempdir(), "test.xlsx")
  
  featureWeightMatrixFile <- system.file("extdata", "Morenos-Mateo.csv", package = "CRISPRseek")
  gRNAs <- findgRNAs(inputFilePath = system.file("extdata", "cpf1.fa", package = "CRISPRseek"), 
                     findPairedgRNAOnly = FALSE,
                     pairOutputFile = pairOutputFile, 
                     PAM = "NGG",
                     PAM.location = "3prime",
                     PAM.size = 3, 
                     overlap.gRNA.positions = c(17, 18), 
                     baseBeforegRNA = 4, 
                     baseAfterPAM = 3, 
                     calculategRNAEfficacy = TRUE,
                     rule.set = "CRISPRscan", 
                     featureWeightMatrixFile = featureWeightMatrixFile,
                     efficacyFile = efficacyFile)
  
  gRNAs_truth <- readRDS(system.file("extdata", "gRNAsSpCas9.RDS", package = "CRISPRseek"))
  expect_equal(as.character(gRNAs_truth), as.character(gRNAs))
    
  effi_truth <- readRDS(system.file("extdata", "CRISPRscanEff.RDS", package = "CRISPRseek"))
  effi <- read.xlsx(efficacyFile)
  
  expect_equal(effi_truth, effi)
})
