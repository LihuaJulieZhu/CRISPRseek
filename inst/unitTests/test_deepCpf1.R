test_deepCpf1 <- function() {
   extendedSequence =
        c('GTTATTTGAGCAATGCCACTTAATAAACATGTAA',
          'TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT',
          'GTTATTTGAGCAATGCCACTTAATAAACATGTAA',
          'TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT')
   chrom_acc = c(0,1, 1, 0)

   seq.DeepScore = c(53.46983, 55.69931, 53.46983, 55.69931)/100
   DeepScore =  c(7.932924, 46.077484, 44.095499,  8.181222)/100

   effi_noacc  <- deepCpf1(extendedSequence)
   effi_acc <- deepCpf1(extendedSequence, chrom_acc)

   checkEqualsNumeric(effi_acc, DeepScore, tolerance = 0.0001)
   checkEqualsNumeric(effi_noacc, seq.DeepScore, tolerance = 0.0001)
   
   extendedSequence <- readRDS(
   system.file("extdata/all.gRNAs.df.RDS", package = "CRISPRseek"))
   extendedSequence <- extendedSequence[, 5]

   chrom_acc <- rep(c(0,1), 10)

   effi_noacc2  <- deepCpf1(extendedSequence)
   effi_acc2 <- deepCpf1(extendedSequence, chrom_acc)

   seq.DeepScore2  <- deepCpf1(extendedSequence[-c(1,2)])
   DeepScore2 <- deepCpf1(extendedSequence[-c(1,2)], chrom_acc[-c(1,2)])

   checkEqualsNumeric(effi_acc2[-c(1,2)], DeepScore2, tolerance = 0.0001)
   checkEqualsNumeric(effi_noacc2[-c(1,2)], seq.DeepScore2, tolerance = 0.0001)
   checkEqualsNumeric(deepCpf1(extendedSequence[1:3]), c(NA, NA, 0.01009442), tolerance = 0.0001)
   checkEqualsNumeric(deepCpf1(extendedSequence[3]), 0.0100944, tolerance = 0.0001)
   suppressWarnings(checkEquals(deepCpf1(extendedSequence[1]), NA))
}

test_deepCpf1BaseMissing <- function() {
   extendedSequence =
        c('GTTATTTGAGCAATGCCACTTAATAAACATGTAA',
          'TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT',
          'GTTATTTGAGCAATGCCACTTAATAAACATGTAA',
          'TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT')
   new.seq <- c(extendedSequence, "GTAGTTTTGAGAGATTTTATTGGTATTGATTTGT")
   eff <- deepCpf1(new.seq) 

   seq.DeepScore = c(53.46983, 55.69931, 53.46983, 55.69931)/100
  
   checkEqualsNumeric(eff[1:4], seq.DeepScore, tolerance = 0.001)

   extendedSeq <- readRDS(system.file("extdata/testSeqForDeepCpf1.RDS",
           package = "CRISPRseek"))

   eff.40 <- deepCpf1(extendedSeq[1:40])
   eff <- deepCpf1(extendedSeq)
   checkEqualsNumeric(eff.40, eff[1:40])
}

library("BSgenome.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

test_somegRNAsNotInGenomeCpf1 <- function()
{
   gRNAs_truth <- readRDS(system.file("extdata",
         "gRNAsCpf1.RDS", package = "CRISPRseek"))

   results1 <- offTargetAnalysis(inputFilePath = gRNAs_truth,
          findgRNAs = FALSE,
          gRNAoutputName = "testCpf1gRNA.fa",
           annotatePaired = FALSE,
	   BSgenomeName = Hsapiens, chromToSearch = "chrX",
           txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
           orgAnn = org.Hs.egSYMBOL, max.mismatch = 0,    
           baseBeforegRNA = 8, baseAfterPAM = 26,
           rule.set = "DeepCpf1", PAM.size = 4,  gRNA.size = 20,
           overlap.gRNA.positions = c(19, 23), 
		outputDir = getwd(), 
		overwrite = TRUE, PAM.location = "5prime",
             PAM = "TTTN", PAM.pattern = "^TNNN", allowed.mismatch.PAM = 2,
             subPAM.position = c(1,2))

    effi_truth <- deepCpf1(results1$summary[!is.na(results1$summary[,2]), 3])
    checkEquals(as.numeric(results1$summary[!is.na(results1$summary[,2]), 4]),
          effi_truth, tolerance = 0.01)

    effi_truth <- deepCpf1(results1$offtarget$extendedSequence)
    checkEquals(results1$offtarget$gRNAefficacy,
        effi_truth, tolerance = 0.01)
}


test_gRNAInGenomeCpf1 <- function()
{
   results2 <- offTargetAnalysis(inputFilePath = system.file("extdata",
         "cpf1-2.fa", package = "CRISPRseek"), findgRNAsWithREcutOnly =  FALSE,
          findPairedgRNAOnly = FALSE,
           annotatePaired = FALSE,
           BSgenomeName = Hsapiens,
	   chromToSearch = "chr8",     
           txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
           orgAnn = org.Hs.egSYMBOL, max.mismatch = 1,
           baseBeforegRNA = 8, baseAfterPAM = 26,
           rule.set = "DeepCpf1",
           overlap.gRNA.positions = c(19, 23), 
		useEfficacyFromInputSeq = FALSE,
		outputDir = getwd(), 
		overwrite = TRUE, PAM.location = "5prime",PAM.size = 4,
             PAM = "TTTN", PAM.pattern = "^TNNN", allowed.mismatch.PAM = 2,
             subPAM.position = c(1,2))
   
   effi_truth <- deepCpf1(results2$summary[!is.na(results2$summary[,2]), 3])

   checkEquals(as.numeric(results2$summary[!is.na(results2$summary[,2]), 4]), 
      effi_truth, tolerance = 0.01)

   effi_truth <- deepCpf1(results2$offtarget$extendedSequence)

   checkEquals(results2$offtarget$gRNAefficacy, 
      effi_truth, tolerance = 0.01)
}

test_gRNAWithOffCpf1 <- function()
{
   results3 <- offTargetAnalysis(inputFilePath = system.file("extdata",
         "cpf1-2.fa", package = "CRISPRseek"), findgRNAsWithREcutOnly =  FALSE,
          findPairedgRNAOnly = FALSE,
           annotatePaired = FALSE,
           BSgenomeName = Hsapiens,
           chromToSearch = paste0("chr", 1:21),   
           txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
           orgAnn = org.Hs.egSYMBOL, max.mismatch = 3,
           baseBeforegRNA = 8, baseAfterPAM = 26,
           rule.set = "DeepCpf1",
           overlap.gRNA.positions = c(19, 23),
                useEfficacyFromInputSeq = FALSE,
                outputDir = getwd(),
                overwrite = TRUE, PAM.location = "5prime",PAM.size = 4,
             PAM = "TTTN", PAM.pattern = "^TNNN", allowed.mismatch.PAM = 2,
             subPAM.position = c(1,2))
  
   effi_truth <- deepCpf1(results3$summary[!is.na(results3$summary[,2]), 3])

   checkEquals(as.numeric(results3$summary[!is.na(results3$summary[,2]), 4]),
      effi_truth, tolerance = 0.01)

   effi_truth <- deepCpf1(results3$offtarget$extendedSequence)

   checkEquals(results3$offtarget$gRNAefficacy,
      effi_truth, tolerance = 0.01)
}

