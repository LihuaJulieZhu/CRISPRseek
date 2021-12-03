#' @importFrom BiocGenerics subset unlist lapply rbind
#' @importFrom hash hash values
#' @importFrom utils read.csv

getOfftargetScore2 <-
    function(featureVectors, 
        subPAM.activity = hash( AA =0,
          AC =   0,
          AG = 0.259259259,
          AT = 0,
          CA = 0,
          CC = 0,
          CG = 0.107142857,
          CT = 0,
          GA = 0.069444444,
          GC = 0.022222222,
          GG = 1,
          GT = 0.016129032,
          TA = 0,
          TC = 0,
          TG = 0.038961039,
          TT = 0),
        mismatch.activity.file = system.file("extdata", 
            "NatureBiot2016SuppTable19DoenchRoot.csv", 
            package = "CRISPRseek")
)
{
    mismatch.activity <- read.csv(mismatch.activity.file)
    required.col <- c("Mismatch.Type", "Position", "Percent.Active")
    if (length(intersect(colnames(mismatch.activity), required.col)) != 
        length(required.col))  
    {
         stop("Please rename the mismatch activity file column to contain at least
              these 3 column names: Mismatch.Type, Position, Percent.Active\n")
    } 
    position <- mismatch.activity$Position
    r.nu.d.nu <- mismatch.activity$Mismatch.Type
    weights <- mismatch.activity$Percent.Active
    #mismatch.activity[mismatch.activity$Mismatch.Type == "rA:dG" & 
    #    mismatch.activity$Position == 10,]$Percent.Active
    ##### by default weights is a column vector
    ##### the mismatch activity  is given as pos 20, 19, 18,....1 distance from PAM,    ##### Position named as 1, 2, 3, .....20 though 
    ##### and the featureVectors is in the same order now
    ##### so no need to reverse any more. weights = rev(weights)
    featureVectors <- subset(featureVectors, !grepl("N", subPAM))
    featureVectors$score <- 1
    featureVectors$score <- as.numeric(values(
        subPAM.activity, keys = as.character(featureVectors$subPAM))) 
    fv.geThan1Mismatch <- subset(featureVectors, as.numeric(as.character(
        featureVectors$n.mismatch)) >= 1)
    fv.lessThan1Mismatch <- subset(featureVectors, as.numeric(as.character(
        featureVectors$n.mismatch)) < 1)
    if (dim(fv.geThan1Mismatch)[1] > 0)
    {
        mismatch.pos <- fv.geThan1Mismatch[,grep("IsMismatch.pos", 
            colnames(fv.geThan1Mismatch))]
        mismatch.pos <- apply(mismatch.pos, 1, as.numeric) 
        pos <- grep("IsMismatch.pos", colnames(fv.geThan1Mismatch))
        min.pos <- min(pos)
        score.new <- unlist(lapply(1:dim(fv.geThan1Mismatch)[1], function(i)
        {
            mismatch.index <- pos[fv.geThan1Mismatch[i,pos] == 1] - min.pos + 1
            thisMismatch <- unlist(strsplit(as.character(
                fv.geThan1Mismatch[i,]$mismatch.type), ","))
            score.new1 <- fv.geThan1Mismatch[i, ]$score
            for (j in 1:length(mismatch.index))
                score.new1 <- score.new1 *
                    weights[position == 
                        mismatch.index[j] & 
                        r.nu.d.nu == thisMismatch[j]]
            score.new1
        }))
        fv.geThan1Mismatch$score <- score.new
    }
    if (dim(fv.geThan1Mismatch)[1] > 0)
    {
        score <- fv.geThan1Mismatch
        if (dim(fv.lessThan1Mismatch)[1] > 0)
            score <- rbind(fv.lessThan1Mismatch, fv.geThan1Mismatch)
    }
    else
    {
        score <- fv.lessThan1Mismatch
    }
    score$alignment <- as.character(score$alignment)
    score$score <- round(score$score, 6)
    score  <- score[order(c(score$name,score$score),decreasing = TRUE), ]
    unique(score[!is.na(score$score), ])
}
