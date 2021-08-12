test_findgRNAsSpCas9 <- function()
{
   featureWeightMatrixFile <- system.file("extdata", "Morenos-Mateo.csv",
                 package = "CRISPRseek")

   gRNAs <- findgRNAs(inputFilePath = system.file("extdata",
         "cpf1.fa", package = "CRISPRseek"),
         findPairedgRNAOnly=FALSE,
         pairOutputFile = "testpairedgRNAs-SpCas9.xls",
         PAM="NGG", PAM.location = "3prime", PAM.size = 3,
         overlap.gRNA.positions = c(17, 18),
         baseBeforegRNA = 4, baseAfterPAM = 3,
         calculategRNAEfficacy= TRUE,
         rule.set = "CRISPRscan",
         featureWeightMatrixFile = featureWeightMatrixFile,
         efficacyFile = "testCRISPRscanEfficacy.xls")

   gRNAs_truth <- readRDS(system.file("extdata",
         "gRNAsSpCas9.RDS", package = "CRISPRseek"))

   checkEquals(gRNAs, gRNAs_truth)

   effi_truth <- readRDS( system.file("extdata",
         "CRISPRscanEff.RDS", package = "CRISPRseek"))

   effi <- read.table("testCRISPRscanEfficacy.xls", sep = "\t", header = TRUE)

   checkEquals(effi, effi_truth)
}

test_findgRNAsCpf1 <- function()
{
   gRNAs <- findgRNAs(inputFilePath = system.file("extdata",
         "cpf1.fa", package = "CRISPRseek"),
         findPairedgRNAOnly=FALSE,
         pairOutputFile = "testpairedgRNAs-cpf1.xls",
         PAM="TTTN", PAM.location = "5prime", PAM.size = 4,
         overlap.gRNA.positions = c(19,23),
         baseBeforegRNA = 8, baseAfterPAM = 26,
         calculategRNAEfficacy= TRUE,
         rule.set = "DeepCpf1",
         efficacyFile = "testcpf1Efficacy.xls")

   gRNAs_truth <- readRDS(system.file("extdata",
         "gRNAsCpf1.RDS", package = "CRISPRseek"))
 
   checkEquals(gRNAs, gRNAs_truth)

   effi_truth <- readRDS( system.file("extdata",
         "Cpf1Eff.RDS", package = "CRISPRseek")) 
  
   effi <- read.table("testcpf1Efficacy.xls", sep = "\t", header = TRUE)
   
   checkEquals(effi, effi_truth)
}
