
test <- 0
if (test == 1)
{
   potential.PEs <- designPEs(inputSeq,
                PAM.size = PAM.size,
		PAM.location = PAM.location,
                gRNA.size = gRNA.size,
                overlap.gRNA.positions = overlap.gRNA.positions,
                PBS.length = PBS.length,
                paired.gRNAs = paired.gRNAs, 
                RT.template.length = RT.template.length,
                RT.template.pattern = RT.template.pattern,
                corrected.seq = corrected.seq,
                targeted.seq.length.change = targeted.seq.length.change,
                bp.after.target.end = bp.after.target.end,
                target.start = target.start,
                target.end = target.end,
                primeEditingPaired.output =  primeEditingPaired.output)
}

designPEs <- function(inputSeq,
        PAM.size = 3L,
        PAM.location = "3prime",
        gRNA.size = 20L,
        overlap.gRNA.positions = c(17L,18L),
        PBS.length = 13L,
        paired.gRNAs,
        RT.template.length,
        RT.template.pattern,
        corrected.seq,
        targeted.seq.length.change,
	bp.after.target.end = 15L,
	target.start,
	target.end,
        primeEditingPaired.output = "pairedgRNAsPE.xls",
        append = FALSE, col.names = TRUE)
{
   ###Hsap_GATA1_ex2_gR7r or Hsap_GATA1_ex2_gR17f
   ###Obtain the cut.start from the names of gRNAs
   ### insert just before target.start
   ### mutate/delete all bases between target.start and target.end inclusive

        cut.start.plus <- as.numeric(unlist(lapply(strsplit(as.character(paired.gRNAs$ForwardgRNAName), "_"), 
             function(x){gsub("f", "", gsub("gR", "", x[length(x)]))})))
        cut.start.minus <- as.numeric(unlist(lapply(strsplit(as.character(paired.gRNAs$ReversegRNAName), "_"),
             function(x){gsub("r", "", gsub("gR", "", x[length(x)]))})))

	gRNAs.plus.cut.d.target.start <- cut.start.plus - target.start
        gRNAs.minus.cut.d.target.end <- target.end - cut.start.minus
              
        actual.RT.template.length.plus <- target.start - cut.start.plus  - 1 +  
              bp.after.target.end 

        actual.RT.template.length.minus <-  cut.start.minus - target.end -1  + 
              bp.after.target.end

        .cutSiteIn5primeOfTargetSite <- function(cutDistanceFromTarget) {
       		unlist(lapply(cutDistanceFromTarget,
                   function(i) {ifelse(i < 0, TRUE, FALSE)})) } 

    #### need to add if PAM.location is 5prime yet #####
        ForwardgRNA.PAM.start <- cut.start.plus + gRNA.size - min(overlap.gRNA.positions) + 1
        ForwardgRNA.PAM.end <-  ForwardgRNA.PAM.start + as.numeric(PAM.size) - 1
        ReversegRNA.PAM.start <- cut.start.minus - gRNA.size + min(overlap.gRNA.positions) - 1
        ReversegRNA.PAM.end <-  ReversegRNA.PAM.start - as.numeric(PAM.size) + 1
        paired.gRNAs <- cbind(paired.gRNAs,
              ReversegRNA.cut.start = cut.start.minus, 
              ForwardgRNA.cut.start = cut.start.plus,	
              ReversegRNA.PAM.start = ReversegRNA.PAM.start,
              ReversegRNA.PAM.end =  ReversegRNA.PAM.end,
              ForwardgRNA.PAM.start = ForwardgRNA.PAM.start,
              ForwardgRNA.PAM.end = ForwardgRNA.PAM.end,
              ReversegRNA.cut.start.d.targetEnd = -gRNAs.minus.cut.d.target.end,
              ForwardgRNA.cut.start.d.targetStart =  gRNAs.plus.cut.d.target.start,
              ReversegRNA.cut.5prime.targetEnd = 
                 .cutSiteIn5primeOfTargetSite(gRNAs.minus.cut.d.target.end),
              ForwardgRNA.cut.5prime.targetStart =
                 .cutSiteIn5primeOfTargetSite(gRNAs.plus.cut.d.target.start),
              ReversegRNA.RT.template.length = actual.RT.template.length.minus,
              ForwardgRNA.RT.template.length = actual.RT.template.length.plus)
 
        paired.gRNAs <- subset(paired.gRNAs, 
            (ReversegRNA.cut.5prime.targetEnd &
               ReversegRNA.RT.template.length  <= max(RT.template.length) & 
               ReversegRNA.RT.template.length >= min(RT.template.length) &
               (ForwardgRNA.PAM.end < target.start | ForwardgRNA.PAM.start > target.end)) |
            (ForwardgRNA.cut.5prime.targetStart &
               ForwardgRNA.RT.template.length <= max(RT.template.length) &
               ForwardgRNA.RT.template.length >= min(RT.template.length) &
               (ReversegRNA.PAM.start < target.start |ReversegRNA.PAM.end > target.end))
            )
 
        PBS.plus <- reverseComplement(.getSeq(inputSeq, 
            starts = paired.gRNAs$ForwardgRNA.cut.start - PBS.length + 1,
            len = rep(PBS.length, length(paired.gRNAs$ForwardgRNA.cut.start)),
            strand = "plus"))
        PBS.minus <- reverse(.getSeq(inputSeq, 
            starts = paired.gRNAs$ReversegRNA.cut.start, 
            len = rep(PBS.length, length(paired.gRNAs$ReversegRNA.cut.start)),
             strand = "plus"))

        actual.RT.template.plus.seq1 <- .getSeq(inputSeq, 
	      starts = paired.gRNAs$ForwardgRNA.cut.start + 1, 
              len = target.start - paired.gRNAs$ForwardgRNA.cut.start - 1,
              strand = "plus")
        
        actual.RT.template.minus.seq1 <- .getSeq(inputSeq, 
              starts = paired.gRNAs$ReversegRNA.cut.start - 1,  
              len = paired.gRNAs$ReversegRNA.cut.start - target.end - 1,
              strand = "minus")

        if (targeted.seq.length.change < 0)
        {
	       actual.RT.template.plus.seq2 = .getSeq(inputSeq, 
                     starts = target.end + 1, 
                     len = bp.after.target.end, "plus") 
               actual.RT.template.minus.seq2 = .getSeq(inputSeq, 
                     starts = target.start - 1, 
                     len = bp.after.target.end, "minus")
	}
        else if (targeted.seq.length.change == 0)
        {
                paired.gRNAs$ReversegRNA.RT.template.length <- 
                    paired.gRNAs$ReversegRNA.RT.template.length + target.end - target.start + 1 
                paired.gRNAs$ForwardgRNA.RT.template.length <-
		     paired.gRNAs$ForwardgRNA.RT.template.length + target.end - target.start + 1
                actual.RT.template.plus.seq2 <- paste(corrected.seq,
			.getSeq(inputSeq, starts = target.end + 1, 
                            len = bp.after.target.end, strand = "plus"),
                          sep = "")

		actual.RT.template.minus.seq2 = paste(corrected.seq,
                        .getSeq(inputSeq, starts = target.start - 1,
                            len = bp.after.target.end, strand = "minus"),
                          sep = "")
        }
        else
        {
                paired.gRNAs$ReversegRNA.RT.template.length <- 
                    paired.gRNAs$ReversegRNA.RT.template.length + targeted.seq.length.change
                paired.gRNAs$ForwardgRNA.RT.template.length <- 
                     paired.gRNAs$ForwardgRNA.RT.template.length + targeted.seq.length.change
		actual.RT.template.plus.seq2 = paste( corrected.seq,
			.getSeq(inputSeq, starts = target.start,  
                            len = bp.after.target.end, strand = "plus"),
                            sep = "")
                actual.RT.template.minus.seq2 = paste(corrected.seq,
                        .getSeq(inputSeq, starts = target.end,
                            len = bp.after.target.end, strand = "minus"),
                            sep = "")
        }
        actual.RT.template.seq.plus <- reverseComplement(DNAStringSet(paste(actual.RT.template.plus.seq1, 
               actual.RT.template.plus.seq2, sep ="")))
        actual.RT.template.seq.minus <- reverse(DNAStringSet(paste(actual.RT.template.minus.seq2, 
               actual.RT.template.minus.seq1, sep ="")))
        paired.gRNAs <- cbind(paired.gRNAs, 
            ReversegRNA.RT.template.seq = actual.RT.template.seq.minus,
            ForwardgRNA.RT.template.seq = actual.RT.template.seq.plus, 
            ReversegRNA.PBS = PBS.minus,
            ForwardgRNA.PBS = PBS.plus)

        write.table(paired.gRNAs, file = primeEditingPaired.output, sep = "\t", 
            row.names = FALSE, col.names = col.names, append = append) 

        #list(paired.gRNAs, actual.RT.template.plus.seq1, actual.RT.template.plus.seq2, actual.RT.template.minus.seq1,
        #    actual.RT.template.minus.seq2) 
       paired.gRNAs
}

.getSeq <- function(inputSeq, starts, len, strand) {
      seq.len <- nchar(as.character(inputSeq))

      n.rec <- length(starts)
      ends <- unlist(lapply(1:n.rec, function(i) {
	      ifelse(strand=="minus",
                   max(1, min(seq.len, starts[i])),
                   max(1, min(seq.len, starts[i] + len[i] - 1)))}))
      starts <-  unlist(lapply(1:n.rec, function(i) {
              ifelse(strand=="minus",
                 max(1, starts[i] - len[i] + 1),
                 max(1, starts[i]))}))
      seq <-  DNAStringSet(Views(inputSeq,
                 start = starts,
                 end = ends))
      seq
      #### output sequence in plus strand disregard whether input is plus or minus
}
