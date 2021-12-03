test_that("test_getOfftargetScore2", {
    featureVectors <- data.frame(
        IsMismatch.pos1 = c(rep(0,5),1, rep(0,3)),
        IsMismatch.pos2 = rep(0,9),
        IsMismatch.pos3 = rep(0,9), 
        IsMismatch.pos4 = c(rep(0,2),1, rep(0,6)),
        IsMismatch.pos5 = rep(0,9),
        IsMismatch.pos6 = c(0, 1, rep(0,7)),
        IsMismatch.pos7 = rep(0,9),
        IsMismatch.pos8 = c(rep(0,3),1,1, rep(0,4)),
        IsMismatch.pos9 = c(rep(0,4),1,1, rep(0,3)),
        IsMismatch.pos10 = rep(0,9),
        IsMismatch.pos11 = rep(0,9),
        IsMismatch.pos12 = c(rep(0,8), 1),
        IsMismatch.pos13 = c(1, rep(0,7), 1),
        IsMismatch.pos14 = c(rep(0,6), 1, rep(0,2)),
        IsMismatch.pos15 = c(rep(0,6), 1, rep(0,2)),
        IsMismatch.pos16 = c(rep(0,7), 1, 0),
        IsMismatch.pos17 = c(rep(0,7), 1, 0),
        IsMismatch.pos18 = rep(0,9),
        IsMismatch.pos19 = rep(0,9),
        IsMismatch.pos20 = rep(0,9),
        n.mismatch = c(1,1,1,1, 2,2,2,2,2),
        mismatch.type = c("rG:dT", "rA:dC", "rC:dA", "rC:dA", "rU:dC,rC:dA",
            "rU:dT,rC:dA", "rG:dT,rA:dG", "rG:dT,rC:dT", "rU:dC,rC:dA"),
       subPAM = rep("GG",9), score = rep(1,9), alignment = rep("", 9))

   correct.scores <- c(0.9230769, 0.7142857, 0.8000000, 0.8750000, 0.6416667,
        0.8750000, 0.2045455, 0.4666667, 0.35) 
   featureVectors$score <- correct.scores
   correct.scores  <- featureVectors[order(
       c(featureVectors$name,featureVectors$score),decreasing = TRUE), ]
   correct.scores <-  unique(correct.scores[!is.na(correct.scores$score), ])
   featureVectors <- correct.scores

    scores <- CRISPRseek:::getOfftargetScore2(featureVectors)
    expect_equal(featureVectors$score, scores$score, tolerance = 0.001)

    fvNinPAM <- readRDS(system.file("extdata", "featureVecorsWithNinPAM.RDS",
            package = "CRISPRseek"))
    cat("checking N in PAM...")
    expect_equal(dim(CRISPRseek:::getOfftargetScore2(fvNinPAM))[1], 1)
})



