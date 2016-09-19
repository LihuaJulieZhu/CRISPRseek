library(CRISPRseek)
inputFile1Path <- system.file("extdata", "rs362331T.fa",
    package = "CRISPRseek")
inputFile2Path <- system.file("extdata", "rs362331C.fa",
    package = "CRISPRseek")
dev.mode <- 0
if (dev.mode == 1)
{
    seqs1CFD <- compare2Sequences.bk(inputFile1Path, inputFile2Path,
        outputDir = getwd(), 
        overwrite = TRUE,scoring.method = "CFDscore")

    seqs1 <- compare2Sequences.bk(inputFile1Path, inputFile2Path,
        outputDir = getwd(), 
        overwrite = TRUE)

    seqs3CFD <- compare2Sequences3(inputFile1Path, inputFile2Path,
        outputDir = getwd(), 
        overwrite = TRUE,scoring.method = "CFDscore")

    seqs3 <- compare2Sequences3(inputFile1Path, inputFile2Path,
        outputDir = getwd(), 
        overwrite = TRUE)

   seqs1.5prime <- compare2Sequences.bk(inputFile1Path, inputFile2Path,
        outputDir = getwd(), PAM.location = "5prime", PAM = "GCT", PAM.pattern = "^NCN",
        overwrite = TRUE, subPAM.position = c(1,2))

   seqs3.5prime<- compare2Sequences3(inputFile1Path, inputFile2Path,
        outputDir = getwd(), PAM.location = "5prime", PAM = "GCT", PAM.pattern = "^NCN",
        overwrite = TRUE, subPAM.position = c(1,2))
}

seqs2CFD.new <- compare2Sequences(inputFile1Path, inputFile2Path,
        outputDir = getwd(), 
        overwrite = TRUE,scoring.method = "CFDscore")

seqs2.new <- compare2Sequences(inputFile1Path, inputFile2Path,
        outputDir = getwd(), 
        overwrite = TRUE)

seqs2.5prime.new <- compare2Sequences(inputFile1Path, inputFile2Path,
        outputDir = getwd(), PAM.location = "5prime", PAM = "GCT", PAM.pattern = "^NCN",
        overwrite = TRUE, subPAM.position = c(1,2))

save(seqs2.5prime.new, file ="~/seqs2.5prime.new.RData")
save(seqs2.new, file = "~/seqs2.new.RData")
save(seqs2CFD.new, file = "~/seqs2CFD.new.RData")

load(system.file("extdata", "seqs2.RData", package = "CRISPRseek"))
load(system.file("extdata", "seqs2CFD.RData", package = "CRISPRseek"))
load(system.file("extdata", "seqs2.5prime.RData", package = "CRISPRseek"))

test_compare2Sequences <- function() {
cat("testing compare2Sequences \n")
if (!checkEquals(seqs2[,-6], seqs2.new[,-6])){
   cat("testing failed for seqs2 and seqs2.new\n!")
}

if (!checkEquals(seqs2CFD[,-6], seqs2CFD.new[,-6]))
{
   cat("testing failed for seqs2CFD and seqs2CFD.new\n!")
}

if (!checkEquals(seqs2.5prime[,-6], seqs2.5prime.new[,-6]))
{
   cat("testing failed for seqs2.5prime and seqs2.5prime.new\n!")
}

if (dev.mode == 1)
{
    if (checkEquals(seqs2, seqs3)){1}

    if (checkEquals(seqs2CFD, seqs3CFD)){1}

    if (checkEquals(seqs2[, -c(2:4)], seqs1[,-c(2:4)])){1}

    if (checkEquals(seqs2CFD[, -c(2:4)], seqs1CFD[,-c(2:4)])){1}

    if (checkEquals(seqs3[, -c(2:4)], seqs1[,-c(2:4)])){1}

    if (checkEquals(seqs3CFD[, -c(2:4)], seqs1CFD[,-c(2:4)])){1}

    if (checkEquals(seqs2.5prime, seqs3.5prime)){1}

    if (checkEquals(seqs3.5prime[, -c(2:4)], seqs1.5prime[,-c(2:4)])){1}
}
}
