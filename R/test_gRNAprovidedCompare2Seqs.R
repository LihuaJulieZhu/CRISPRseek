library(CRISPRseek)
inputFile2Path <- system.file("extdata", "rs362331C.fa",
    package = "CRISPRseek")

seqs2 <- readRDS(system.file("extdata", "seqs2.RDS",
    package = "CRISPRseek"))

inputgRNA <- DNAStringSet("GTAGATGAGGGAGCAGGCGT")

gRNAs1Provided <- compare2Sequences(inputgRNA, inputFile2Path,
        outputDir = getwd(), findgRNAs = c(FALSE, FALSE), searchDirection = "1to2",
         overwrite = TRUE)


gRNAs2Provided <- compare2Sequences( inputFile2Path, inputgRNA,
        outputDir = getwd(), findgRNAs = c(FALSE, FALSE), searchDirection = "2to1",
         overwrite = TRUE)

test_gRNAprovidedCompare2Seqs <- function() {
   for (i in 2:11)
   {
       if (!checkEquals(as.character(gRNAs1Provided[1, i]), as.character(seqs2[1,i])))
         cat(paste(i, "th column gRNAs1 provided differ from the first row of seqs2"))
   }
    if (!checkEquals(gRNAs1Provided[, -c(3,4,7,8,14)], gRNAs2Provided[, -c(3,4,7,8,14)]))
    {
        cat("gRNAs1 provided differ from gRNAs2 provided")
    }

    if (!checkEquals(as.character(gRNAs1Provided[,3]), as.character(gRNAs2Provided[,4])))
    { 
        cat("3rd column gRNAs1 provided differ from 4th column in gRNAs2 provided")
    }

    if (!checkEquals(as.character(gRNAs1Provided[,4]), as.character(gRNAs2Provided[,3])))
    { 
        cat("4th column gRNAs1 provided differ from 3rd column in gRNAs2 provided")
    }

    if (!checkEquals(as.character(gRNAs1Provided[,7]), as.character(gRNAs2Provided[,8])))
    { 
       cat("7th column gRNAs1 provided differ from 8th column in gRNAs2 provided")
    }

    if (!checkEquals(as.character(gRNAs1Provided[,8]), as.character(gRNAs2Provided[,7])))
    { 
       cat("8th column gRNAs1 provided differ from 7th column in gRNAs2 provided")
    }

    if (!checkEquals(abs(as.numeric(as.character(gRNAs1Provided[,14]))), abs(as.numeric(as.character(gRNAs2Provided[,14])))))
    { 
       cat("the absolute value of 14th column gRNAs1 provided differ from 14th column in gRNAs2 provided")
    }
}
