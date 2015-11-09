.getgRNA.cut.sites <- function(subject, subjectname, PAM ="NGG", 
    gRNA.pattern = "", gRNA.size = 20, cut.site = 17,
    PAM.size = 3, calculategRNAEfficacy = TRUE, baseBeforegRNA = 4,
    baseAfterPAM = 3,
    reverse.subject = FALSE)
{
     seq.len <- nchar(as.character(subject))
     pos.PAMs <- unlist(gregexpr(PAM, subject, perl = TRUE,
         ignore.case = TRUE, fixed = FALSE))
     pos.PAMs <- pos.PAMs[pos.PAMs != -1]
     starts.gRNA <- pos.PAMs - gRNA.size
     ends.gRNA <- pos.PAMs - 1
     pos.PAMs <- pos.PAMs[starts.gRNA > 0]
     ends.gRNA <- ends.gRNA[starts.gRNA > 0]
     starts.gRNA <- starts.gRNA[starts.gRNA > 0] 
    
     if (gRNA.pattern != "")
     {
         gRNA.seqs <- as.character(Views(subject, 
             start = starts.gRNA,
             end = ends.gRNA))
         pos.set2 <- unlist(gregexpr(gRNA.pattern, gRNA.seqs,
             perl = TRUE, ignore.case = TRUE, fixed = FALSE))
         pos.PAMs <- pos.PAMs[pos.set2 == 1]
     }
     seq <- as.character(Views(subject,
         start = starts.gRNA,
         end = ends.gRNA + PAM.size))
     extendedSequence <- seq
     if (calculategRNAEfficacy)
     {
         extended.starts <- starts.gRNA - baseBeforegRNA
         extended.starts[extended.starts < 0] <- 1
         extended.ends <- ends.gRNA + PAM.size + baseAfterPAM
         extended.ends[extended.ends > length(subject)] <- length(subject)
         extendedSequence <- as.character(Views(subject, start = extended.starts,
             end = extended.ends))
     }
     if (reverse.subject)
     {
         gRNAs.cut <- cbind(seq, paste( subjectname,"_gR",
             (seq.len - (starts.gRNA + cut.site -1 ) + 1), "r", sep = ""),
             (seq.len - starts.gRNA + 1), "-", extendedSequence)
     }
     else
     {
         gRNAs.cut <-
             cbind(seq, paste(subjectname,"_gF", (starts.gRNA + cut.site - 1),
             "f", sep = ""), starts.gRNA, "+", extendedSequence)
     }
     gRNAs.cut
}
findgRNAs <-
    function (inputFilePath, format = "fasta", PAM = "NGG", PAM.size = 3,
        findPairedgRNAOnly = FALSE, annotatePaired = TRUE,
        enable.multicore = FALSE, n.cores.max = 6,
        gRNA.pattern = "", gRNA.size = 20, 
	overlap.gRNA.positions = c(17,18),
        min.gap = 0, max.gap = 20, pairOutputFile, name.prefix = "",
	    featureWeightMatrixFile = system.file("extdata", 
        "DoenchNBT2014.csv", package = "CRISPRseek"), baseBeforegRNA = 4, 
	    baseAfterPAM = 3,
	    calculategRNAEfficacy = FALSE, efficacyFile) 
{
    if (missing(inputFilePath)) {
        stop("inputFilePath containing the searching sequence 
	   or a DNAStringSet object is required!")
    }
    cut.site <- min(overlap.gRNA.positions)
    if (class(inputFilePath) != "DNAStringSet")
    {
    	if (! file.exists(inputFilePath)) {
           stop("inputfile specified as ", inputFilePath, " does not exists!")
        }
        if (format == "fasta" || format == "fastq")
        {
            subjects <- readDNAStringSet(inputFilePath, format, use.names = TRUE)
        }
        else
        {
            stop("format needs to be either fasta or fastq !")
        }
    }
    else
    {
	subjects <- inputFilePath
    }
    PAM <- translatePattern(PAM)
    PAM <- paste("(?=", PAM, ")", sep="")

    gRNA.pattern <- translatePattern(gRNA.pattern)
    min.subject <- gRNA.size + PAM.size 
    subjects <-  subjects[width(subjects) >= min.subject,]
    if (length(subjects) == 0)
    {
        stop("The input file contains no sequence! This could be caused by 
            wrong format of the file. If file is created in mac, you could 
            reformat to text by typing tr \"\\r\" \"\\n\" >newfile in the 
            command line")
    }
    toAppend = FALSE
    colNames = TRUE
    names(subjects) <- gsub("'", "", names(subjects))
    names(subjects) <- gsub(" ", "", names(subjects))
    names(subjects) <- gsub("\t", ":", names(subjects))
     
    all.gRNAs.df <- do.call(rbind, lapply(1:length(subjects), function(p)
   {
        subject <- subjects[[p]]
        subjectname <- names(subjects)[p]
        revsubject <- reverseComplement(subject)
        plus.gRNAs <- 
           .getgRNA.cut.sites(subject, subjectname, PAM = PAM, 
               gRNA.pattern = gRNA.pattern, 
               gRNA.size = gRNA.size,
               cut.site = cut.site,
               PAM.size = PAM.size, 
               calculategRNAEfficacy = calculategRNAEfficacy,
               baseBeforegRNA = baseBeforegRNA,
               baseAfterPAM = baseAfterPAM)
       
       minus.gRNAs <-
           .getgRNA.cut.sites(revsubject, subjectname, PAM = PAM,
               gRNA.pattern = gRNA.pattern,
               gRNA.size = gRNA.size,
               cut.site = cut.site,
               PAM.size = PAM.size,
               calculategRNAEfficacy = calculategRNAEfficacy,
               baseBeforegRNA = baseBeforegRNA,
               baseAfterPAM = baseAfterPAM,
               reverse.subject = TRUE) 

          plus.index <- numeric()
          minus.index <- numeric()
          forward.index <- 0
          reverse.index <- 0
          if (length(plus.gRNAs) > 1)
            n.plus.gRNAs <- dim(plus.gRNAs)[1]
          else
            n.plus.gRNAs <- 0
          if (length(minus.gRNAs) > 1)
            n.minus.gRNAs <- dim(minus.gRNAs)[1]
          else
            n.minus.gRNAs <- 0
        if (findPairedgRNAOnly && (n.minus.gRNAs  * n.plus.gRNAs) == 0)
        {
            warning(paste("No paired gRNAs found in the input sequence", subjectname))
            all.gRNAs <- DNAStringSet()
        }
	else if (annotatePaired || findPairedgRNAOnly)
	{
            temp = matrix(nrow =0, ncol=5)
            colnames(temp)[1:5] <- c( "ReversegRNAPlusPAM",
                "ReversegRNAName", "ForwardgRNAPlusPAM",
                "ForwardgRNAName", "gap")
	        if (p >1)
	        {
	            toAppend = TRUE
	            colNames = FALSE
	        }
            write.table(temp, file = pairOutputFile, sep = "\t",
                row.names = FALSE, quote = FALSE, append = toAppend, 
	            col.names = colNames)
            if (n.minus.gRNAs > 0 && n.plus.gRNAs > 0)
            {
                for (j in 1:n.plus.gRNAs)
                {
                    for (k in 1:n.minus.gRNAs)
                    {
                        if ((as.numeric(as.character(plus.gRNAs[j,3])) - 1 - 
                            as.numeric(as.character(minus.gRNAs[k,3]))) >= min.gap
                            && (as.numeric(as.character(plus.gRNAs[j,3])) - 
                            as.numeric(as.character(minus.gRNAs[k,3]))) <= 
                            max.gap)
                        {
                            plus.index <- c(plus.index, j)
                            minus.index <- c(minus.index, k)
                        } ### if paired
                    } #for k
                }## for j
                if (length(minus.index) > 0 && length(plus.index) > 0)
                {
                    paired <- cbind(minus.gRNAs[minus.index,1], 
                        minus.gRNAs[minus.index, 2], plus.gRNAs[plus.index,1], 
                        plus.gRNAs[plus.index, 2],
                        gap = as.numeric(as.character(
                        plus.gRNAs[plus.index,3])) - 1 - as.numeric(
                        as.character(minus.gRNAs[minus.index, 3])))
                    colnames(paired)[1:5] <- c( "ReversegRNAPlusPAM", 
                        "ReversegRNAName", "ForwardgRNAPlusPAM", 
                        "ForwardgRNAName", "gap")
                    if(findPairedgRNAOnly)
                    {
                        plus.index <- unique(plus.index)
                        minus.index <- unique(minus.index)
                        all.gRNAs <- DNAStringSet(c(plus.gRNAs[plus.index, 1], 
                            minus.gRNAs[minus.index,1]))
                        names(all.gRNAs) <- c(plus.gRNAs[plus.index,2], 
                            minus.gRNAs[minus.index,2])
                    }
                    if (dim(paired)[1] == 1)
                        write.table(paired, file = pairOutputFile, sep = "\t", 
                            row.names = FALSE, quote = FALSE, append = toAppend,
			                col.names = colNames)
                    else
                        write.table(paired[order(as.character(paired[,4])), ], 
                            file = pairOutputFile, sep = "\t", row.names = FALSE,
                            quote = FALSE, append = toAppend, col.names = colNames)
                }
                else if (findPairedgRNAOnly)
                {
                    warning(paste("No paired gRNAs found for sequence",
                        subjectname))
                   all.gRNAs <- DNAStringSet()
                }
            }### if plus.gRNAs and minus.gRNAs not empty
            if (! findPairedgRNAOnly)
            {
                all.gRNAs <- DNAStringSet(c(plus.gRNAs[,1], minus.gRNAs[,1]))
                if (length(all.gRNAs) == 0)
                    warning(paste("No gRNAs found in the input sequence", 
                        subjectname))	
            }
	 } ### annotatePaired and (no paired found or findPairedOnly)   
         else
         {
             all.gRNAs <- DNAStringSet(c(plus.gRNAs[,1], minus.gRNAs[,1]))
             if (length(all.gRNAs) == 0)
                 warning(paste("No gRNAs found in the input sequence",
                     subjectname))
         }
	 forEffi <- rbind(plus.gRNAs, minus.gRNAs)
	 forEffi <- subset(forEffi, forEffi[,1] %in% as.character(all.gRNAs))	
	 if (length(all.gRNAs) >0)
	     forEffi
    })) ## do call subjects
    if (calculategRNAEfficacy)
    {
        featureWeightMatrix <- read.csv(featureWeightMatrixFile, header=TRUE)
        effi <- calculategRNAEfficiency(all.gRNAs.df[,5], 
	    baseBeforegRNA = baseBeforegRNA, 
	    featureWeightMatrix = featureWeightMatrix, 
            enable.multicore = enable.multicore,
            n.cores.max = n.cores.max,
            gRNA.size = gRNA.size) 
        extendedSequences <- cbind(all.gRNAs.df, effi)
        colnames(extendedSequences)  <- c("gRNAplusPAM", "name", "start", "strand", 
	    "extendedSequence", "gRNAefficacy")
        extendedSequences[nchar(extendedSequences[,5]) 
	    < baseBeforegRNA + gRNA.size + PAM.size + baseAfterPAM, 6] <- 
                "extended sequence too short"
	write.table(extendedSequences,
             file = efficacyFile, sep="\t", row.names = FALSE)
    }
    #else
    #{
#	extendedSequences <- all.gRNAs.df[,1:4]
#	colnames(extendedSequences)  = c("gRNAplusPAM", "name", "start", "strand")
#    }
    all.gRNAs <- DNAStringSet(all.gRNAs.df[,1])
    names(all.gRNAs) = all.gRNAs.df[,2]
    #list(all.gRNAs=all.gRNAs, extendedSequences = extendedSequences)
    all.gRNAs
}
