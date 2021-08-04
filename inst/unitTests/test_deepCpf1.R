 in_df <-  data.frame('gRNA' = paste0("gRNA", 1:4), 
    sequence =
        c('GTTATTTGAGCAATGCCACTTAATAAACATGTAA',
          'TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT',
          'GTTATTTGAGCAATGCCACTTAATAAACATGTAA',
          'TGACTTTGAATGGAGTCGTGAGCGCAAGAACGCT'),
        chrom = c(0,1, 1, 0))

out_truth <- cbind(in_df, 
    seq.DeepScore = c(53.46983, 55.69931, 53.46983, 55.69931),
    DeepScore =  c(7.932924, 46.077484, 44.095499,  8.181222))

out_df  <- deepCpf1(in_df)

test_deepCpf1 <- function() {
   checkEqualsNumeric(out_df$DeepScore, out_truth$DeepScore, tolerance = 0.0001)
   checkEqualsNumeric(out_df$seq.DeepScore, out_truth$seq.DeepScore, tolerance = 0.0001)
   checkEquals(out_df[,1:3], out_truth[,1:3])
}
