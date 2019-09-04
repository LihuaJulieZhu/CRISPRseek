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
    allGenes <- genes(txdb, columns = "gene_id", single.strand.genes.only=TRUE)
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
