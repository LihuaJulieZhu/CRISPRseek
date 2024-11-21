# This is a legit unit test only if the seqs2, seqs2CFD, seqs2.5prime are golden truth.

test_that("test_compare2Sequences", {
  inputFile1Path <- system.file("extdata", "rs362331T.fa", package = "CRISPRseek")
  inputFile2Path <- system.file("extdata", "rs362331C.fa", package = "CRISPRseek")
  seqs2CFD.new <- compare2Sequences(inputFile1Path, 
                                    inputFile2Path, 
                                    outputDir = tempdir(), 
                                    overwrite = TRUE,
                                    scoring.method = "CFDscore")
  seqs2.new <- compare2Sequences(inputFile1Path, 
                                 inputFile2Path,
                                 outputDir = tempdir(), 
                                 overwrite = TRUE)
  seqs2.5prime.new <- compare2Sequences(inputFile1Path,
                                        inputFile2Path,
                                        outputDir = tempdir(), 
                                        PAM.location = "5prime", 
                                        PAM = "GCT", 
                                        PAM.pattern = "^NCN",
                                        overwrite = TRUE, 
                                        subPAM.position = c(1,2))
  
  seqs2 <- readRDS(system.file("extdata", 
                               "seqs2.RDS", 
                               package = "CRISPRseek"))
  seqs2CFD <- readRDS(system.file("extdata", 
                                  "seqs2CFD.RDS", 
                                  package = "CRISPRseek"))
  seqs2.5prime <- readRDS(system.file("extdata", 
                                      "seqs2.5prime.RDS", 
                                      package = "CRISPRseek"))
  
  expect_equal(seqs2[, -6], seqs2.new[order(as.numeric(rownames(seqs2.new))), -6], tolerance = 0.001)
  expect_equal(seqs2CFD[, -6], seqs2CFD.new[order(as.numeric(rownames(seqs2CFD.new))), -6], tolerance = 0.001)
  expect_equal(seqs2.5prime[, -6], seqs2.5prime[order(as.numeric(rownames(seqs2.5prime.new))), -6], tolerance = 0.001)
})

