test_that(
   "extendedSequence is required for DeepCpf1 algorithm",
   {
      expect_error(deepCpf1())
   })

test_that(
   "deepcpf1 expects the sequence to be 34 base long",
   {
      expect_warning(deepCpf1("AATGAAA"))
   }
)

test_that("test_deepCpf1", {
    extendedSequence = c("GTTATTTGAGCAATGCCACTTAATAAACATGTAA", 
        "TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT", "GTTATTTGAGCAATGCCACTTAATAAACATGTAA", 
        "TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT")
    chrom_acc = c(0, 1, 1, 0)
    seq.DeepScore = c(53.46983, 55.69931, 53.46983, 55.69931)/100
    DeepScore = c(7.932924, 46.077484, 44.095499, 8.181222)/100
    effi_noacc <- deepCpf1(extendedSequence)
    effi_acc <- deepCpf1(extendedSequence, chrom_acc)
    expect_equal(DeepScore, effi_acc, tolerance = 0.001)
    expect_equal(seq.DeepScore, effi_noacc, tolerance = 0.001)
    extendedSequence <- readRDS(system.file("extdata/all.gRNAs.df.RDS", 
        package = "CRISPRseek"))
    extendedSequence <- extendedSequence[, 5]
    chrom_acc <- rep(c(0, 1), 10)
    effi_noacc2 <- deepCpf1(extendedSequence)
    effi_acc2 <- deepCpf1(extendedSequence, chrom_acc)
    seq.DeepScore2 <- deepCpf1(extendedSequence[-c(1, 2)])
    DeepScore2 <- deepCpf1(extendedSequence[-c(1, 2)], chrom_acc[-c(1, 
        2)])
    expect_equal(DeepScore2, effi_acc2[-c(1, 2)], tolerance = 0.001)
    expect_equal(seq.DeepScore2, effi_noacc2[-c(1, 2)], tolerance = 0.001)
    expect_equal(c(NA, NA, 0.01009442), deepCpf1(extendedSequence[1:3]), tolerance = 0.001)
    expect_equal(0.0100944, deepCpf1(extendedSequence[3]), tolerance = 0.0001)
    suppressWarnings(expect_equal(deepCpf1(extendedSequence[1]), 
        NA))
})

test_that("test_deepCpf1BaseMissing", {
    extendedSequence = c("GTTATTTGAGCAATGCCACTTAATAAACATGTAA", 
        "TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT", "GTTATTTGAGCAATGCCACTTAATAAACATGTAA", 
        "TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT")
    new.seq <- c(extendedSequence, "GTAGTTTTGAGAGATTTTATTGGTATTGATTTGT")
    eff <- deepCpf1(new.seq)
    seq.DeepScore = c(53.46983, 55.69931, 53.46983, 55.69931)/100
    expect_equal(seq.DeepScore, eff[1:4], tolerance = 0.001)
    extendedSeq <- readRDS(system.file("extdata/testSeqForDeepCpf1.RDS", 
        package = "CRISPRseek"))
    eff.40 <- deepCpf1(extendedSeq[1:40])
    eff <- deepCpf1(extendedSeq)
    expect_equal(eff[1:40], eff.40, tolerance = 0.001)
})

test_that("test_gRNAInGenomeCpf1", {
    results2 <- offTargetAnalysis(inputFilePath = system.file("extdata", 
        "cpf1-2.fa", package = "CRISPRseek"), findgRNAsWithREcutOnly = FALSE, 
        findPairedgRNAOnly = FALSE, annotatePaired = FALSE, BSgenomeName = Hsapiens, 
        chromToSearch = "chr8", txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
        orgAnn = org.Hs.egSYMBOL, max.mismatch = 1, baseBeforegRNA = 8, 
        baseAfterPAM = 26, rule.set = "DeepCpf1", overlap.gRNA.positions = c(19, 
            23), useEfficacyFromInputSeq = FALSE, outputDir = getwd(), 
        overwrite = TRUE, PAM.location = "5prime", PAM.size = 4, 
        PAM = "TTTN", PAM.pattern = "^TNNN", allowed.mismatch.PAM = 2, 
        subPAM.position = c(1, 2))
    effi_truth <- deepCpf1(results2$summary[!is.na(results2$summary[, 
        2]), 3])
    expect_equal(effi_truth, as.numeric(results2$summary[!is.na(results2$summary[, 
        2]), 4]), tolerance = 0.001)
    effi_truth <- deepCpf1(results2$offtarget$extendedSequence)
    expect_equal(effi_truth, results2$offtarget$gRNAefficacy, tolerance = 0.001)
})

test_that("test_gRNAWithOffCpf1", {
    results3 <- offTargetAnalysis(inputFilePath = system.file("extdata", 
        "cpf1-2.fa", package = "CRISPRseek"), findgRNAsWithREcutOnly = FALSE, 
        findPairedgRNAOnly = FALSE, annotatePaired = FALSE, BSgenomeName = Hsapiens, 
        chromToSearch = paste0("chr", 1:21), txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
        orgAnn = org.Hs.egSYMBOL, max.mismatch = 3, baseBeforegRNA = 8, 
        baseAfterPAM = 26, rule.set = "DeepCpf1", overlap.gRNA.positions = c(19, 
            23), useEfficacyFromInputSeq = FALSE, outputDir = getwd(), 
        overwrite = TRUE, PAM.location = "5prime", PAM.size = 4, 
        PAM = "TTTN", PAM.pattern = "^TNNN", allowed.mismatch.PAM = 2, 
        subPAM.position = c(1, 2))
    effi_truth <- deepCpf1(results3$summary[!is.na(results3$summary[, 
        2]), 3])
    expect_equal(effi_truth, as.numeric(results3$summary[!is.na(results3$summary[, 
        2]), 4]), tolerance = 0.001)
    effi_truth <- deepCpf1(results3$offtarget$extendedSequence)
    expect_equal(effi_truth, results3$offtarget$gRNAefficacy, tolerance = 0.001)
})

test_that("test_somegRNAsNotInGenomeCpf1", {
    gRNAs_truth <- readRDS(system.file("extdata", "gRNAsCpf1.RDS", 
        package = "CRISPRseek"))
    results1 <- offTargetAnalysis(inputFilePath = gRNAs_truth, 
        findgRNAs = FALSE, gRNAoutputName = "testCpf1gRNA.fa", 
        annotatePaired = FALSE, BSgenomeName = Hsapiens, chromToSearch = "chrX", 
        txdb = TxDb.Hsapiens.UCSC.hg19.knownGene, orgAnn = org.Hs.egSYMBOL, 
        max.mismatch = 0, baseBeforegRNA = 8, baseAfterPAM = 26, 
        rule.set = "DeepCpf1", PAM.size = 4, gRNA.size = 20, 
        overlap.gRNA.positions = c(19, 23), outputDir = getwd(), 
        overwrite = TRUE, PAM.location = "5prime", PAM = "TTTN", 
        PAM.pattern = "^TNNN", allowed.mismatch.PAM = 2, subPAM.position = c(1, 
            2))
    effi_truth <- deepCpf1(results1$summary[!is.na(results1$summary[, 
        2]), 3])
    expect_equal(effi_truth, as.numeric(results1$summary[!is.na(results1$summary[, 
        2]), 4]), tolerance = 0.001)
    effi_truth <- deepCpf1(results1$offtarget$extendedSequence)
    expect_equal(effi_truth, results1$offtarget$gRNAefficacy, tolerance = 0.001)
})
