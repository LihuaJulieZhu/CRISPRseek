#library(BSgenome.Hsapiens.UCSC.hg19)
#library(CRISPRseek)
### T2test100bp.bed
#offTargetAnalysisOfPeakRegions(gRNA = "T2.fa",
#   peaks = "T2test100bp.bed",
#   format=c("fasta", "bed"),
#   peaks.withHeader = FALSE, BSgenomeName = Hsapiens,
#   upstream = 50, downstream =50, PAM.size = 3, gRNA.size = 20,
#   PAM = "NGG", PAM.pattern = "NNN$", max.mismatch = 8,
#   outputDir="TS2bedplusminusMerged",
#   allowed.mismatch.PAM = 3, overwrite = TRUE
#   )

offTargetAnalysisOfPeakRegions <- function(gRNA, peaks, 
   format=c("fasta", "bed"),
   peaks.withHeader = FALSE, BSgenomeName,
   upstream = 50, downstream =50, PAM.size = 3, gRNA.size = 20,
   PAM = "NGG", PAM.pattern = "NNN$", max.mismatch = 8,
   outputDir, allowed.mismatch.PAM = 3, overwrite = TRUE
   )
{
   thePeaks <- read.table(peaks, sep="\t", header = peaks.withHeader, 
      stringsAsFactors = FALSE)
   if (format[2] == "bed")
   {
     if (dim(thePeaks)[2] >= 4 )
       colnames(thePeaks)[1:4] <- c("chromosome", "peak_start", "peak_end", "names")
     if(dim(thePeaks)[2] >= 5)
       colnames(thePeaks)[5] = "peak_score"
     if(dim(thePeaks)[2] >= 6)
       colnames(thePeaks)[6] = "peak_strand"
   }
   else
   {
      stop("only bed file with at least 4 columns are supported")
   }
   thePeaks[, 4] <- gsub(" ", "", thePeaks[,4])
   TS2 <- compare2Sequences(inputFile1Path = gRNA, inputFile2Path = peaks, 
      findgRNAsWithREcutOnly = FALSE, format = format, header = peaks.withHeader,
      BSgenomeName = BSgenomeName, upstream = upstream, downstream = downstream,
      minREpatternSize = 4, searchDirection = "1to2",
      overlap.gRNA.positions = c(17, 18), findPairedgRNAOnly = FALSE, 
      min.gap = 0, max.gap = 20, gRNA.name.prefix = "gRNA", PAM.size = PAM.size, 
      gRNA.size = 20, PAM = PAM, PAM.pattern = PAM.pattern, max.mismatch = max.mismatch,
      outputDir = outputDir, foldgRNAs = FALSE,
      allowed.mismatch.PAM = allowed.mismatch.PAM, overwrite = overwrite
     )
   TS2 <- cbind(names = as.character(unlist(lapply(as.character(TS2$offTarget), 
      function(temp) {paste(strsplit(temp,":")[[1]][1:6], collapse=":")}))), TS2)
   excluding.columns = which(colnames(TS2) %in% 
      c("scoreForSeq1", "targetInSeq1", "gRNAefficacy", "scoreDiff"))
   TS2 <- TS2[, -excluding.columns]
   colnames(TS2)[colnames(TS2) == "scoreForSeq2"] = "predicted_cleavage_score"
   colnames(TS2)[colnames(TS2) == "targetInSeq2"] = "offTarget_sequence"
   offTargetOffset <- do.call(rbind, lapply(as.character(TS2$offTarget),
      function(temp) 
      { 
         as.numeric(strsplit(strsplit(temp,":")[[1]][7], "-")[[1]][1:2])
         }))
   offtargets <- merge(TS2, thePeaks, by = "names")
   offTargetOffset[, 1] <- offTargetOffset[, 1] - upstream + offtargets$peak_start
   offTargetOffset[, 2] <- offTargetOffset[, 2] - upstream + offtargets$peak_start
   offtargets <- cbind(offtargets, offTarget_Start = offTargetOffset[1], 
      offTarget_end = offTargetOffset[2])
   offtargets <- merge(offtargets, thePeaks, all = TRUE)
   
   write.table(offtargets, file = "offTargetsInPeakRegions.xls", sep="\t", row.names = FALSE)
   list(offtargets = offtargets, TS2 = TS2, thePeaks = thePeaks)
}
