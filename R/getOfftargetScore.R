#' Calculate score for each off target
#' 
#' Calculate score for each off target with given feature vectors and weights
#' vector
#' 
#' score is calculated using the weights and algorithm by Hsu et al., 2013
#' cited in the reference section
#' 
#' @param featureVectors a data frame generated from
#' buildFeatureVectorForScoring. It contains 
#' \itemize{
#' \item{IsMismatch.posX} - {Indicator variable
#' indicating whether this position X is mismatch or not, (1 means yes and 0
#' means not). X takes on values from 1 to gRNA.size, representing all positions in the guide RNA (gRNA).}
#' \item{strand} - {strand of the off target, + for plus and - for minus strand}
#' \item{chrom} - {chromosome of the off target}
#' \item{chromStart} - {start position of the off
#' target}
#' \item{chromEnd} - {end position of the off target}
#' \item{name} - {gRNA name}
#' \item{gRNAPlusPAM} - {gRNA sequence with PAM sequence concatenated}
#' \item{OffTargetSequence} - {the genomic sequence of the off target}
#' \item{n.mismatch} - {number of mismatches between the off target and the gRNA}
#' \item{forViewInUCSC} - {string for viewing in UCSC genome browser, e.g., chr14:31665685-31665707}
#' \item{score} - {score of the off target}
#' \item{mismatche.distance2PAM} - {a comma separated
#' distances of all mismatches to PAM, e.g., 14,11 means one mismatch is 14 bp
#' away from PAM and the other mismatch is 11 bp away from PAM} 
#' \item{alignment} - {alignment between gRNA and off target, e.g., ......G..C.......... means
#' that this off target aligns with gRNA except that G and C are mismatches}
#' \item{NGG} - {this off target contains canonical PAM or not, 1 for yes and 0 for no}
#' \item{mean.neighbor.distance.mismatch} - {mean distance between neighboring
#' mismatches}
#' }
#' @param weights a numeric vector size of gRNA length, default c(0, 0, 0.014,
#' 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732,
#' 0.828, 0.615, 0.804, 0.685, 0.583) which is used in Hsu et al., 2013 cited
#' in the reference section
#' @return a data frame containing 
#' \itemize{
#' {strand} - {strand of the match, + for plus and
#' - for minus strand}
#' \item{chrom} - {chromosome of the off target}
#' \item{chromStart} - {start
#' position of the off target}
#' \item{chromEnd} - {end position of the off target}
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
#' that this off target aligns with gRNA except that G and C are mismatches}
#' \item{NGG} - {this off target contains canonical PAM or not, 1 for yes and 0 for no}
#' \item{mean.neighbor.distance.mismatch} - {mean distance between neighboring
#' mismatches}
#' }
#' @note %% ~~further notes~~
#' @author Lihua Julie Zhu
#' @seealso offTargetAnalysis
#' @references Patrick D Hsu, David A Scott, Joshua A Weinstein, F Ann Ran,
#' Silvana Konermann, Vineeta Agarwala, Yinqing Li, Eli J Fine, Xuebing Wu,
#' Ophir Shalem, Thomas J Cradick, Luciano A Marraffini, Gang Bao & Feng Zhang
#' (2013) DNA targeting specificity of rNA-guided Cas9 nucleases. Nature
#' Biotechnology 31:827-834
#' @keywords misc
#' @examples
#' 
#'     hitsFile <-  system.file("extdata", "hits.txt", 
#'         package = "CRISPRseek")
#'     hits <- read.table(hitsFile, sep = "\t", header = TRUE,
#'         stringsAsFactors = FALSE)
#'     featureVectors <- buildFeatureVectorForScoring(hits)
#'     getOfftargetScore(featureVectors)
#' @importFrom BiocGenerics subset unlist lapply rbind
#' @export
getOfftargetScore <-
    function(featureVectors, 
        weights=c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 
        0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583)) 
{
    ##### by default weights is a column vector
    ##### the weights is given as pos 20, 19, 18,....1 distance from PAM, 
    ##### and the featureVectors is in the same order now
    ##### so no need to reverse any more. weights = rev(weights)
    fv.lessThan2Mismatch <- subset(featureVectors, as.numeric(as.character(
        featureVectors$n.mismatch)) < 2)
    if (dim(fv.lessThan2Mismatch)[1] > 0)
    {
        mismatch.pos <- fv.lessThan2Mismatch[, grep("IsMismatch.pos", 
            colnames(fv.lessThan2Mismatch))]
        mismatch.pos = apply(mismatch.pos, 1, as.numeric) 
        fv.lessThan2Mismatch$score =  as.numeric(100 * 
            (1 - weights %*% mismatch.pos))
    }
    fv.geThan2Mismatch <- subset(featureVectors, as.numeric(as.character(
        featureVectors$n.mismatch)) >= 2)
    if (dim(fv.geThan2Mismatch)[1] > 0)
    {
        mismatch.pos <- fv.geThan2Mismatch[,grep("IsMismatch.pos", 
            colnames(fv.geThan2Mismatch))]
        mismatch.pos <- apply(mismatch.pos, 1, as.numeric) 
        fv.geThan2Mismatch$score <-  100
        pos <- grep("IsMismatch.pos", colnames(fv.lessThan2Mismatch))
        min.pos <- min(pos)
        score.new = unlist(lapply(1:dim(fv.geThan2Mismatch)[1], function(i)
        {
            mismatch.index <- pos[fv.geThan2Mismatch[i,pos] == 1] - min.pos + 1
            score.new <- fv.geThan2Mismatch[i, ]$score
            for (j in mismatch.index)
                score.new <- score.new * (1 - weights[j])
            score.new
        }))
        fv.geThan2Mismatch$score <- score.new
        ### Penalize consecutive mismatches more than mismatches far apart
        fv.geThan2Mismatch$score <- fv.geThan2Mismatch$score / 
            (((19-as.numeric(as.character(
            fv.geThan2Mismatch$mean.neighbor.distance.mismatch))) / 19) * 4 + 1)
        ### Penalize even more with more mismatches (quadratic)
        fv.geThan2Mismatch$score <- fv.geThan2Mismatch$score / as.numeric(
            as.character(fv.geThan2Mismatch$n.mismatch)) ^ 2;
        fv.geThan2Mismatch$score[fv.geThan2Mismatch$score < 0] <- 1;
    }
    if (dim(fv.geThan2Mismatch)[1] > 0)
    {
        score <- fv.geThan2Mismatch
        if (dim(fv.lessThan2Mismatch)[1] > 0)
            score <- rbind(fv.lessThan2Mismatch, fv.geThan2Mismatch)
    }
    else
    {
        score <- fv.lessThan2Mismatch
    }
    score$alignment <- as.character(score$alignment)
    score$score <- round(score$score, 1)
    score  <- score[order(c(score$name,score$score),decreasing = TRUE), ]
    unique(score[!is.na(score$score), ])
}
