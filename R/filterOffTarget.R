#' filter off targets and generate reports.
#' 
#' filter off targets that meet the criteria set by users such as minimum
#' score, topN. In addition, off target was annotated with flank sequence, gRNA
#' cleavage efficiency and whether it is inside an exon or not if fetchSequence
#' is set to TRUE and annotateExon is set to TRUE
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param scores a data frame output from getOfftargetScore. It contains 
#' \itemize{
#' \item{strand} - {strand of the off target, + for plus and - for minus strand}
#' \item{chrom} - {chromosome of the off target}
#' \item{chromStart} - {start position of the offtarget}
#' \item{chromEnd} - {end position of the offtarget}
#' \item{name} - {gRNA name} 
#' \item{gRNAPlusPAM} - {gRNA sequence with PAM sequence concatenated}
#' \item{OffTargetSequence} - {the genomic sequence of the off target}
#' \item{n.mismatch} - {number of mismatches between the off target and the gRNA}
#' \item{forViewInUCSC} - {string for viewing in UCSC genome browser, e.g., chr14:31665685-31665707}
#' \item{score} - {score of the off target}
#' \item{mismatch.distance2PAM} - {a comma separated
#' distances of all mismatches to PAM, e.g., 14,11 means one mismatch is 14 bp
#' away from PAM and the other mismatch is 11 bp away from PAM}
#' \item{alignment} - {alignment between gRNA and off target, e.g., ......G..C.......... means
#' that this off target aligns with gRNA except that G and C are
#' mismatches}
#' \item{NGG} - {this off target contains canonical PAM or not, 1 for yes
#' and 0 for no)} 
#' \item{mean.neighbor.distance.mismatch} - {mean distance between
#' neighboring mismatches}
#' }
#' @param min.score minimum score of an off target to included in the final
#' output, default 0.5
#' @param topN top N off targets to be included in the final output, default
#' 100
#' @param topN.OfftargetTotalScore top N off target used to calculate the total
#' off target score, default 10
#' @param annotateExon Choose whether or not to indicate whether the off target
#' is inside an exon or not, default TRUE
#' @param txdb TxDb object, for creating and using TxDb object, please refer to
#' GenomicFeatures package. For a list of existing TxDb object, please search
#' for annotation package starting with Txdb at
#' http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData,
#' such as TxDb.Rnorvegicus.UCSC.rn5.refGene for rat,
#' TxDb.Mmusculus.UCSC.mm10.knownGene for mouse,
#' TxDb.Hsapiens.UCSC.hg19.knownGene for human,
#' TxDb.Dmelanogaster.UCSC.dm3.ensGene for Drosophila and
#' TxDb.Celegans.UCSC.ce6.ensGene for C.elegans
#' @param orgAnn organism annotation mapping such as org.Hs.egSYMBOL in
#' org.Hs.eg.db package for human
#' @param ignore.strand default to TRUE
#' @param outputDir the directory where the off target analysis and reports
#' will be written to
#' @param oneFilePergRNA write to one file for each gRNA or not, default to
#' FALSE
#' @param fetchSequence Fetch flank sequence of off target or not, default TRUE
#' @param upstream upstream offset from the off target start, default 200
#' @param downstream downstream offset from the off target end, default 200
#' @param BSgenomeName BSgenome object. Please refer to available.genomes in
#' BSgenome package. For example, 
#' \itemize{
#' \item{BSgenome.Hsapiens.UCSC.hg19} - {for hg19}
#' \item{BSgenome.Mmusculus.UCSC.mm10} - {for mm10}
#' \item{BSgenome.Celegans.UCSC.ce6} - {for ce6}
#' \item{BSgenome.Rnorvegicus.UCSC.rn5} - {for rn5}
#' \item{BSgenome.Dmelanogaster.UCSC.dm3} - {for dm3}
#' }
#' @param baseBeforegRNA Number of bases before gRNA used for calculating gRNA
#' efficiency, default 4
#' @param baseAfterPAM Number of bases after PAM used for calculating gRNA
#' efficiency, default 3
#' @param gRNA.size The size of the gRNA, default 20 for spCas9
#' @param PAM.location PAM location relative to gRNA. For example, spCas9 PAM
#' is located on the 3 prime while cpf1 PAM is located on the 5 prime
#' @param PAM.size PAM length, default 3 for spCas9
#' @param featureWeightMatrixFile Feature weight matrix file used for
#' calculating gRNA efficiency. By default DoenchNBT2014 weight matrix is used.
#' To use alternative weight matrix file, please input a csv file with first
#' column containing significant features and the second column containing the
#' corresponding weights for the features. Please see Doench et al., 2014 for
#' details.
#' @param rule.set Specify a rule set scoring system for calculating gRNA
#' efficacy.
#' @param chrom_acc Optional binary variable indicating chromatin accessibility 
#' information with 1 indicating accessible and 0 not accessible.
#' @param calculategRNAefficacyForOfftargets Default to TRUE to output gRNA
#' efficacy for offtargets as well as ontargets. Set it to FALSE if only need
#' gRNA efficacy calculated for ontargets only to speed up the analysis. Please
#' refer to https://support.bioconductor.org/p/133538/#133661 for potential use
#' cases of offtarget efficacies.
#' @return \item{offtargets }{a data frame with off target analysis results}
#' \item{summary }{a data frame with summary of the off target analysis
#' results}
#' @note %% ~~further notes~~
#' @author Lihua Julie Zhu
#' @seealso offTargetAnalysis
#' @references Doench JG, Hartenian E, Graham DB, Tothova Z, Hegde M, Smith I,
#' Sullender M, Ebert BL, Xavier RJ, Root DE. Rational design of highly active
#' sgRNAs for CRISPR-Cas9-mediated gene inactivation. Nat Biotechnol. 2014 Sep
#' 3. doi: 10.1038 nbt.3026 Lihua Julie Zhu, Benjamin R. Holmes, Neil Aronin
#' and Michael Brodsky. CRISPRseek: a Bioconductor package to identify
#' target-specific guide RNAs for CRISPR-Cas9 genome-editing systems. Plos One
#' Sept 23rd 2014
#' @keywords misc
#' @examples
#' 
#'     library(CRISPRseek)
#'     library("BSgenome.Hsapiens.UCSC.hg19")
#'     library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'     library(org.Hs.eg.db)
#'     hitsFile <-  system.file("extdata", "hits.txt", package="CRISPRseek")
#'     hits <- read.table(hitsFile, sep = "\t", header = TRUE, 
#'         stringsAsFactors = FALSE)
#'     featureVectors <- buildFeatureVectorForScoring(hits)
#'     scores <- getOfftargetScore(featureVectors)
#'     outputDir <- getwd() 
#'     results <- filterOffTarget(scores, BSgenomeName = Hsapiens, 
#'         txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'          orgAnn = org.Hs.egSYMBOL, outputDir = outputDir, 
#'         min.score = 0.1, topN = 10, topN.OfftargetTotalScore = 5)
#'     results$offtargets
#'     results$summary
#' @importFrom utils read.csv write.table read.table
#' @importFrom BiocGenerics unlist cbind 
#' @importFrom BSgenome getSeq
#' @importFrom S4Vectors merge
#' @importFrom GenomeInfoDb seqlengths
#' @export
filterOffTarget <-
    function(scores, min.score = 0.01, topN = 200, topN.OfftargetTotalScore = 20,
    annotateExon = TRUE, txdb, orgAnn, ignore.strand = TRUE, outputDir, oneFilePergRNA = FALSE,
    fetchSequence = TRUE, upstream = 200, downstream = 200, BSgenomeName,
    baseBeforegRNA = 4, baseAfterPAM = 3, gRNA.size = 20, PAM.location = "3prime", PAM.size = 3,
    featureWeightMatrixFile = system.file("extdata", "DoenchNBT2014.csv", 
	package = "CRISPRseek"),
    rule.set = c("Root_RuleSet1_2014", "Root_RuleSet2_2016", "CRISPRscan", "DeepCpf1"),
    chrom_acc,
   calculategRNAefficacyForOfftargets = TRUE
)
{
    rule.set <- match.arg(rule.set)
   
    if (featureWeightMatrixFile != ""  && file.exists(featureWeightMatrixFile))
	{
		featureWeightMatrix <- read.csv(featureWeightMatrixFile, header=TRUE)
	}
	if (fetchSequence && (missing(BSgenomeName) || 
        class(BSgenomeName) != "BSgenome")) 
    {
        stop("To fetch sequences, BSgenomeName is required as BSgenome object!")
    }
    if (annotateExon && ( missing(txdb) || (class(txdb) != "TxDb" && class(txdb) != "TranscriptDb")))
    {
        stop("To indicate whether an offtarget is inside an exon, txdb is
            required as TxDb object!")
    }
    if (annotateExon && missing(orgAnn))
    {
	    warning("orgAnn was not included. See the updated manual for information about how to use the orgAnn parameter to generate gene identifiers in the offTarget output file")
    }
    scores <- scores[scores$score >= min.score,]
    if (length(grep("IsMismatch.pos", colnames(scores))) > 0)
        scores <- scores[,-c(grep("IsMismatch.pos", colnames(scores)))]
	if (substr(outputDir, nchar(outputDir), nchar(outputDir)) != .Platform$file.sep)
    {
        outputDir <- paste(outputDir, "", sep = .Platform$file.sep)
    }	
    if ( ! file.exists(outputDir))
    {
        dir.create(outputDir)
    }
    OfftargetFile <-paste(outputDir, "OfftargetAnalysis.xls", sep = "")
    OfftargetSummary <-paste(outputDir, "Summary.xls", sep = "")
    gRNAsPlusPAM<- unique(scores$name)
    names <- gRNAsPlusPAM
    top5OfftargetTotalScore <- numeric(length(names))
    topNOfftargetTotalScore <- top5OfftargetTotalScore
    temp <- cbind(names, gRNAsPlusPAM, top5OfftargetTotalScore, 
        topNOfftargetTotalScore)
    mismatch.distance2PAM <- matrix(ncol = 11, nrow = length(names))
    append <- FALSE
    for (i in 1:length(gRNAsPlusPAM))
    {
        this.score <- scores[scores$name == gRNAsPlusPAM[i],]
        this.score <- this.score[order(this.score$score, this.score$n.mismatch, decreasing = c(TRUE, FALSE)),]
        maxN <- min(topN+1, dim(this.score)[1])
        this.score <- this.score[1:maxN,]
        maxN.totalScore <- min(maxN, (topN.OfftargetTotalScore + 1))
        if (this.score$n.mismatch[1] == 0 && as.numeric(as.character(this.score$NGG[1])) == 1)
        {
           start.ind <- 2
           end.ind <- min(maxN, 6)
           end.forSummary <- 11
        }
        else
        {
           start.ind <- 1
           maxN <- maxN - 1
           maxN.totalScore <- maxN.totalScore - 1
           end.forSummary <- 10
           end.ind <- min(maxN, 5) 
        }
        temp[i,3] <- sum(this.score$score[start.ind:end.ind])
        #if (maxN < 6)
        #    temp[i,3] <- sum(this.score$score[2:maxN])
        #else
        #    temp[i,3] <- sum(this.score$score[2:6]) 
        
        if (maxN < maxN.totalScore)
            temp[i,4] <- sum(this.score$score[start.ind:maxN])
	else
            temp[i,4] <- sum(this.score$score[start.ind:maxN.totalScore])
        temp[i,2] <- unique(this.score$gRNAPlusPAM)
        mismatch.distance2PAM[i,] <- 
	   ifelse(as.character(this.score$mismatch.distance2PAM[1]) == "", "NMM",
            "perfect match not found")
        # end.forSummary is 10 if no on-target found, otherwise 11

        forSummary <- this.score[start.ind:end.forSummary,]
	forSummary <- forSummary[order(forSummary$score, decreasing=TRUE),]
	mismatch.distance2PAM[i,2:11] <-
           as.character(forSummary$mismatch.distance2PAM)

	if (dim(forSummary)[1] < 10)
	   mismatch.distance2PAM[i, (dim(forSummary)[1] +1):11] <- "NA"
	this.score <- cbind(name = this.score$name, 
                gRNAPlusPAM = this.score$gRNAPlusPAM,
                OffTargetSequence = this.score$OffTargetSequence,
                score = this.score$score, n.mismatch = this.score$n.mismatch, 
                mismatch.distance2PAM = 
                as.character(this.score$mismatch.distance2PAM), 
                alignment = this.score$alignment,
                NGG = as.character(this.score$NGG),
                forViewInUCSC = this.score$forViewInUCSC, 
                strand = this.score$strand, 
                chrom = this.score$chrom, chromStart = this.score$chromStart,
                chromEnd = this.score$chromEnd)
        if (oneFilePergRNA & dim(this.score)[1] > 0)
            write.table(this.score[!is.na(this.score[,grep("score", 
                colnames(this.score))]),], 
                file = paste( outputDir, "OfftargetAnalysis-", 
                as.character(temp[i,1]), ".xls", sep = ""), sep = "\t",
                row.names = FALSE)
        if (i == 1 && dim(this.score)[1] > 0)
        {
            write.table(this.score, file = OfftargetFile, sep = "\t", 
                row.names = FALSE,  append = append)
            append <- TRUE
        }
        else if (dim(this.score)[1] > 0)
        {
            #this.score <- this.score[!is.na(this.score[,grep("score", 
            ##colnames(this.score))]),]
            write.table(this.score, file = OfftargetFile, sep = "\t",
                row.names = FALSE, col.names = FALSE, append = append)
            append <- TRUE
        }
    }
    temp <- cbind(temp, mismatch.distance2PAM)
    colnames(temp)[5] <- "top1Hit(onTarget)MMdistance2PAM"
    colnames(temp)[4] <- paste("top", topN.OfftargetTotalScore, 
       "OfftargetTotalScore", sep = "")
    colnames(temp)[6:15] <- paste("topOfftarget", 1:10, "MMdistance2PAM",
        sep = "")
    Offtargets <- read.table(OfftargetFile, sep = "\t", header = TRUE, 
        stringsAsFactors = FALSE)
	if (annotateExon)
	{
		Offtargets <- annotateOffTargets(Offtargets, txdb, orgAnn, ignore.strand)
	}
    ontargets <- subset(Offtargets, Offtargets$n.mismatch == 0)
    if (!calculategRNAefficacyForOfftargets && dim(ontargets)[1] > 0)
    {
	chr <- as.character(ontargets$chrom)
        strand <- as.character(ontargets$strand)
        Start <- ifelse(strand=="-",
              as.numeric(as.character( ontargets$chromStart)) - baseAfterPAM,
              as.numeric(as.character( ontargets$chromStart)) - baseBeforegRNA)
        End <- ifelse(strand=="-",
              as.numeric(as.character( ontargets$chromEnd)) + as.numeric(baseBeforegRNA),
              as.numeric(as.character( ontargets$chromEnd)) + as.numeric(baseAfterPAM))
    }
    else if (calculategRNAefficacyForOfftargets && dim(Offtargets)[1] > 0)
    {
	chr <- as.character(Offtargets$chrom)
        strand <- as.character(Offtargets$strand)
        if (PAM.location == "3prime") 
        {
           Start <- ifelse(strand=="-",
              as.numeric(as.character( Offtargets$chromStart)) - baseAfterPAM,
              as.numeric(as.character( Offtargets$chromStart)) - baseBeforegRNA)
           End <- ifelse(strand=="-",
              as.numeric(as.character( Offtargets$chromEnd)) + as.numeric(baseBeforegRNA),
              as.numeric(as.character( Offtargets$chromEnd)) + as.numeric(baseAfterPAM))
        }
        else
        {
	   Start <- ifelse(strand=="-",
              as.numeric(as.character( Offtargets$chromStart)) - baseAfterPAM + gRNA.size,
              as.numeric(as.character( Offtargets$chromStart)) - baseBeforegRNA + PAM.size)
           End <- ifelse(strand=="-",
              as.numeric(as.character( Offtargets$chromEnd)) + as.numeric(baseBeforegRNA) - PAM.size,
              as.numeric(as.character( Offtargets$chromEnd)) + as.numeric(baseAfterPAM) - gRNA.size)
        }
    }
    if ((calculategRNAefficacyForOfftargets && dim(Offtargets)[1] > 0) || (!calculategRNAefficacyForOfftargets && dim(ontargets)[1] > 0))
    {
        starts <- unlist(apply(cbind(Start,1), 1, max))
        ends <- unlist(apply(cbind(End, seqlengths(BSgenomeName)[chr]), 1,min))
        
        extendedSequence <- getSeq(BSgenomeName, names = chr, start = starts,
           end = ends, strand = strand, width = NA, as.character = TRUE)
        if (rule.set == "Root_RuleSet1_2014")
	{
            gRNAefficiency <- calculategRNAEfficiency(extendedSequence, 
		baseBeforegRNA = baseBeforegRNA,
		featureWeightMatrix = featureWeightMatrix)
        }
        else if (rule.set == "Root_RuleSet2_2016")
        {
            gRNAefficiency <- calculategRNAEfficiency2(extendedSequence)
        }
         else if (rule.set == "CRISPRscan")
        {
            gRNAefficiency <- calculategRNAEfficiencyCRISPRscan(extendedSequence,
              featureWeightMatrix = featureWeightMatrix)
        }
        else if (rule.set == "DeepCpf1")
        {
	    gRNAefficiency <- round(deepCpf1(extendedSequence = extendedSequence, 
                 chrom_acc = chrom_acc), 3)
        }
       if (!calculategRNAefficacyForOfftargets && dim(ontargets)[1] > 0)
       {
	   ontargets <- cbind(ontargets,  extendedSequence = extendedSequence, gRNAefficacy = gRNAefficiency)
           Offtargets <- merge(Offtargets, ontargets, all = TRUE)
        }
        else {
           Offtargets  <- cbind(Offtargets,  extendedSequence = extendedSequence, gRNAefficacy = gRNAefficiency)
        }
      }
      if (fetchSequence)
	{
           strand <- as.character(Offtargets$strand)
           chr <- as.character(Offtargets$chrom)
	   Start <- ifelse(strand=="-", 
	      as.numeric(as.character(Offtargets$chromStart)) - as.numeric(downstream),
              as.numeric(as.character(Offtargets$chromStart)) - as.numeric(upstream))
	   End <- ifelse(strand=="-",             
              as.numeric(as.character(Offtargets$chromEnd)) + as.numeric(upstream),
              as.numeric(as.character(Offtargets$chromEnd)) + as.numeric(downstream)) 
	   starts <- unlist(apply(cbind(Start,1), 1, max))
	   ends <- unlist(apply(cbind(End, seqlengths(BSgenomeName)[chr]), 1,min))
	   seq <- getSeq(BSgenomeName, names = chr, start = starts,
                end = ends, strand = strand, width = NA, as.character = TRUE)
	   Offtargets <- cbind(Offtargets, flankSequence = seq)
	}
    colnames(Offtargets)[colnames(Offtargets) == "NGG"] = "isCanonicalPAM"
    write.table(temp, file = OfftargetSummary, sep = "\t", row.names = FALSE)
    write.table(Offtargets[order(as.character(Offtargets$name), 
        -as.numeric(as.character(Offtargets$score)), 
        as.character(Offtargets$OffTargetSequence)),], 
        file = OfftargetFile, sep = "\t", row.names = FALSE)
    list(offtargets = unique(Offtargets), summary = unique(temp))
}
