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
