#' annotate off targets 
#' 
#' Annotate Off targets to indicate whether each one (respectively) is inside an exon or intron, as well as
#' the gene ID if inside the gene.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param scores a data frame output from getOfftargetScore or filterOfftarget.
#' It contains
#' \itemize{
#' \item{strand} - {strand of the off target ((+) for plus and (-) for minus
#' strand)}
#' \item{chrom} - {chromosome of the off target}
#' \item{chromStart} - {start position of
#' the off target}
#'\item{chromEnd} - {end position of the off target}
#'\item{name} - {gRNA
#' name}
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
#' and 0 for no}
#' \item{mean.neighbor.distance.mismatch} - {mean distance between
#' neighboring mismatches}
#' }
#' @param txdb TxDb object. For creating and using TxDb object, please refer to
#' GenomicFeatures package. \\
#' For a list of existing TxDb object, please search
#' for annotation package starting with Txdb at
#' http://www.bioconductor.org/packages/release/BiocViews.html#___AnnotationData,
#' such as 
#' \itemize{
#' \item{TxDb.Rnorvegicus.UCSC.rn5.refGene} - {for rat}
#' \item{TxDb.Mmusculus.UCSC.mm10.knownGene} - {for mouse}
#' \item{TxDb.Hsapiens.UCSC.hg19.knownGene} - {for human}
#' \item{TxDb.Dmelanogaster.UCSC.dm3.ensGene} - {for Drosophila}
#' \item{TxDb.Celegans.UCSC.ce6.ensGene} - {for C.elegans}
#' }
#' @param orgAnn organism annotation mapping such as org.Hs.egSYMBOL. Which lives in the
#' org.Hs.eg.db package for humans.
#' @param ignore.strand default to TRUE
#' 
#' @return a Data Frame with Off Target annotation
#' @note %% ~~further notes~~
#' @author Lihua Julie Zhu
#' @seealso offTargetAnalysis
#' @references Lihua Julie Zhu, Benjamin R. Holmes, Neil Aronin and Michael
#' Brodsky. CRISPRseek: a Bioconductor package to identify target-specific
#' guide RNAs for CRISPR-Cas9 genome-editing systems. Plos One Sept 23rd 2014
#' @keywords misc
#' @examples
#' 
#'     library(CRISPRseek)
#'     #library("BSgenome.Hsapiens.UCSC.hg19")
#'     library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'     library(org.Hs.eg.db)
#'     hitsFile <-  system.file("extdata", "hits.txt", package="CRISPRseek")
#'     hits <- read.table(hitsFile, sep = "\t", header = TRUE, 
#'         stringsAsFactors = FALSE)
#'     featureVectors <- buildFeatureVectorForScoring(hits)
#'     scores <- getOfftargetScore(featureVectors)
#'     outputDir <- getwd() 
#'     results <- annotateOffTargets(scores, 
#'         txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#'          orgAnn = org.Hs.egSYMBOL)
#'     results
#' @importFrom GenomeInfoDb seqlevels 
#' @importFrom BSgenome seqnames 
#' @importFrom IRanges IRanges overlapsAny findOverlaps 
#' @importFrom GenomicRanges GRanges GRangesList 
#' @importFrom BiocGenerics cbind rbind unlist toTable 
#' @importFrom S4Vectors queryHits subjectHits Rle
#' @importFrom methods as
#' @export 
annotateOffTargets <- function(scores, txdb, orgAnn, ignore.strand = TRUE)
{
    score.RD <- GRanges(seqnames = Rle(scores$chrom), 
        ranges = IRanges(start = scores$chromStart, 
       	end = scores$chromEnd, names = scores$forViewInUCSC),
	strand = scores$strand)
    allExons <- as(exons(txdb, columns="gene_id"),"GRanges")
    if (length(grep("Chr",seqnames(allExons))) == 0 && 
		length(grep("Chr", scores$chrom)) >0 )
    {
	seqlevels(allExons) = paste("Chr", seqlevels(allExons), sep="")
    }
    allExons <- allExons[as.character(seqnames(allExons)) %in% 
      	unique(as.character(seqnames(score.RD))),]
    if (ignore.strand)
    {
        score.plus.RD <- GRanges(seqnames = Rle(scores$chrom),
            ranges = IRanges(start = scores$chromStart,
            end = scores$chromEnd, names = scores$forViewInUCSC),
            strand = "+")
        score.minus.RD <- GRanges(seqnames = Rle(scores$chrom),
            ranges = IRanges(start = scores$chromStart,
            end = scores$chromEnd, names = scores$forViewInUCSC),
            strand = "-")
        score.list <- GRangesList(score.plus.RD, score.minus.RD)
        score.RD <- c(score.list[[1]], score.list[[2]])
        scores <- rbind(scores, scores)
     } 
    ann.scores <- overlapsAny(score.RD, allExons, minoverlap = 1L, 
        type = "any",ignore.strand=FALSE)
    inExon <- cbind(forViewInUCSC = names(score.RD),
        inExon = unlist(ann.scores))
    inExon[inExon[,2] == FALSE, 2] <- ""
    allGenes <- suppressMessages(genes(txdb, columns = "gene_id", single.strand.genes.only=TRUE))
    if (length(grep("Chr",seqnames(allGenes))) == 0 && 
                length(grep("Chr", scores$chrom)) >0 )
    {
        seqlevels(allGenes) = paste("Chr", seqlevels(allGenes), sep="")
    }
    overlapGenes <- findOverlaps(score.RD, allGenes, minoverlap = 1L, 
	type = "any",ignore.strand=FALSE)
	query.ind <- queryHits(overlapGenes)
        entrez_id <- character(dim(scores)[1])
        symbol <- entrez_id
	entrez_id[query.ind] <- 
		unlist(allGenes[subjectHits(overlapGenes),]$gene_id)
	entrez_id <- as.character(entrez_id)
	entrez_id[is.na(entrez_id)] = ""
	scores <- cbind(scores, entrez_id = entrez_id, symbol = entrez_id)
	
	if (length(queryHits(overlapGenes)) > 0 && !missing(orgAnn) && 
		class(orgAnn) == "AnnDbBimap" )
	{
	    egSYMBOL <- toTable(orgAnn)
	    if(length(grep("flybase_id", colnames(egSYMBOL)[2])) >0)
            {
                m <- match(scores$entrez_id, egSYMBOL$flybase_id)
		scores$symbol <- egSYMBOL[,1][m]
	    }
	    else
	    {
		m <- match(scores$entrez_id, egSYMBOL$gene_id)
		scores$symbol <- egSYMBOL$symbol[m]
	    }
		scores$symbol[is.na(scores$symbol)] = ""
	}
        scores <- cbind(scores,inExon)
	inIntron <- entrez_id
	inIntron[scores$entrez_id != "" & scores$inExon == ""] = TRUE
	inIntron[scores$entrez_id == "" | scores$inExon == TRUE] = ""
	scores <- cbind(scores, inIntron = inIntron)
        scores <- scores[order(scores$entrez_id, na.last = TRUE, decreasing=TRUE),]
        scores <- scores[!duplicated(scores$forViewInUCSC) | scores$entrez_id != "",]
    if (length(grep("FBgn",entrez_id[1])) > 0 && !missing(orgAnn) && 
		class(orgAnn) == "AnnDbBimap")
    {
	    temp.id <- scores$entrez_id
	    scores$entrez_id <- scores$symbol
	    scores$symbol <- temp.id
	    rm(temp.id)
    }
    if (!missing(orgAnn) && class(orgAnn) == "AnnDbBimap")
    {
        scores <- cbind(name = as.character(scores$name),
    	    gRNAPlusPAM = as.character(scores$gRNAPlusPAM),
            OffTargetSequence = as.character(scores$OffTargetSequence),
            inExon = as.character(scores$inExon),
	        inIntron = as.character(scores$inIntron),
	        entrez_id = as.character(scores$entrez_id),
	        gene = as.character(scores$symbol),
            score = scores$score, n.mismatch = scores$n.mismatch, 
			mismatch.distance2PAM = as.character(
               	scores$mismatch.distance2PAM), 
            alignment = as.character(scores$alignment),
            NGG = as.character(scores$NGG),
            forViewInUCSC = as.character(scores$forViewInUCSC), 
            strand = as.character(scores$strand),
            chrom = as.character(scores$chrom), chromStart = scores$chromStart,
            chromEnd = scores$chromEnd)
    }
    else
    {
        scores <- cbind(name = as.character(scores$name),
			gRNAPlusPAM = as.character(scores$gRNAPlusPAM),
			OffTargetSequence = as.character(scores$OffTargetSequence),
			inExon = as.character(scores$inExon),
			inIntron = as.character(scores$inIntron),
			entrez_id = as.character(scores$entrez_id),
			score = scores$score, n.mismatch = scores$n.mismatch,
			mismatch.distance2PAM = as.character(
				scores$mismatch.distance2PAM), 
			alignment = as.character(scores$alignment),
			NGG = as.character(scores$NGG),
			forViewInUCSC = as.character(scores$forViewInUCSC), 
			strand = as.character(scores$strand),
			chrom = as.character(scores$chrom), chromStart = scores$chromStart,
			chromEnd = scores$chromEnd)
    }
    unique(data.frame(scores))
}
