extendedSequence =
        c('GTTATTTGAGCAATGCCACTTAATAAACATGTAA',
          'TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT',
          'GTTATTTGAGCAATGCCACTTAATAAACATGTAA',
          'TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT')
chrom_acc = c(0,1, 1, 0)

seq.DeepScore = c(53.46983, 55.69931, 53.46983, 55.69931)
DeepScore =  c(7.932924, 46.077484, 44.095499,  8.181222)

effi_noacc  <- deepCpf1(extendedSequence)
effi_acc <- deepCpf1(extendedSequence, chrom_acc)

extendedSequence <- readRDS(
   system.file("extdata/all.gRNAs.df.RDS", package = "CRISPRseek"))
extendedSequence <- extendedSequence[, 5]

chrom_acc <- rep(c(0,1), 10)

effi_noacc2  <- deepCpf1(extendedSequence)
effi_acc2 <- deepCpf1(extendedSequence, chrom_acc)

seq.DeepScore2  <- deepCpf1(extendedSequence[-c(1,2)])
DeepScore2 <- deepCpf1(extendedSequence[-c(1,2)], chrom_acc[-c(1,2)])

test_deepCpf1 <- function() {
   checkEqualsNumeric(effi_acc, DeepScore, tolerance = 0.0001)
   checkEqualsNumeric(effi_noacc, seq.DeepScore, tolerance = 0.0001)
   checkEqualsNumeric(effi_acc2[-c(1,2)], DeepScore2, tolerance = 0.0001)
   checkEqualsNumeric(effi_noacc2[-c(1,2)], seq.DeepScore2, tolerance = 0.0001)
   checkEqualsNumeric(deepCpf1(extendedSequence[1:3]), c(NA, NA, 1.009442), tolerance = 0.0001)
   checkEqualsNumeric(deepCpf1(extendedSequence[3]), 1.00944, tolerance = 0.0001)
   suppressWarnings(checkEquals(deepCpf1(extendedSequence[1]), NA))
}

library("BSgenome.Hsapiens.UCSC.hg19")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

test_offtargetCpf1 <- function()
{
   results1 <- offTargetAnalysis(inputFilePath = system.file("extdata",
         "cpf1.fa", package = "CRISPRseek"), findgRNAsWithREcutOnly =  FALSE,
          findPairedgRNAOnly = FALSE,
           annotatePaired = FALSE,
	   BSgenomeName = Hsapiens, chromToSearch = "chrX",
           txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
           orgAnn = org.Hs.egSYMBOL, max.mismatch = 0,    
           baseBeforegRNA = 8, baseAfterPAM = 26,
           rule.set = "DeepCpf1", PAM.size = 4,  gRNA.size = 20,
           overlap.gRNA.positions = c(19, 23), 
		outputDir = "~/Documents/chrX", 
		overwrite = TRUE, PAM.location = "5prime",
             PAM = "TTTN", PAM.pattern = "^TNNN", allowed.mismatch.PAM = 2,
             subPAM.position = c(1,2))
    

   results2 <- offTargetAnalysis(inputFilePath = system.file("extdata",
         "cpf1.fa", package = "CRISPRseek"), findgRNAsWithREcutOnly =  FALSE,
          findPairedgRNAOnly = FALSE,
           annotatePaired = FALSE,
	      chromToSearch = "",     
           baseBeforegRNA = 8, baseAfterPAM = 26,
           rule.set = "DeepCpf1",
           overlap.gRNA.positions = c(19, 23), 
		useEfficacyFromInputSeq = TRUE,
		outputDir = "~/Documents/", 
		overwrite = TRUE, PAM.location = "5prime",PAM.size = 4,
             PAM = "TTTN", PAM.pattern = "^TNNN", allowed.mismatch.PAM = 2,
             subPAM.position = c(1,2))
}
