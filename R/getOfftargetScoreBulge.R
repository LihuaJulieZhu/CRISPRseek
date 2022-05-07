.updateScore <- function(fv.geThan1, fv.lessThan1,
                         col.prefix = "IsInsertion.pos",
                         weights, position, type)
{
    if (dim(fv.geThan1)[1] > 0)
    {
        indel.pos <- fv.geThan1[,grep(col.prefix,
                                      colnames(fv.geThan1))]
        indel.pos <- apply(indel.pos, 1, as.numeric)
        pos <- grep(col.prefix, colnames(fv.geThan1))
        min.pos <- min(pos)

        if (length(pos) > 0)
        {
            score.new <- unlist(lapply(1:nrow(fv.geThan1), function(i)
            {
                indel.index <- pos[fv.geThan1[i, pos] == 1] - min.pos + 1
                if (col.prefix == "IsInsertion.pos")
                    this.type <- unlist(strsplit(as.character(
                        fv.geThan1[i,]$gRNA.insertion), ","))
                else
                    this.type <- unlist(strsplit(as.character(
                        fv.geThan1[i,]$gRNA.deletion), ","))
                score.new1 <- fv.geThan1[i, ]$score
                for (j in 1:length(indel.index))
                    score.new1 <- score.new1 *
                            weights[position ==
                                indel.index[j] &
                                type == this.type[j]]
                score.new1
            }))
            fv.geThan1$score <- score.new
        }
    }
    if (dim(fv.geThan1)[1] > 0)
    {
        score <- fv.geThan1
        if (dim(fv.lessThan1)[1] > 0)
            score <- rbind(fv.lessThan1, fv.geThan1)
    }
    else
    {
        score <- fv.lessThan1
    }
    score
}

#' @importFrom BiocGenerics subset unlist lapply rbind
#' @importFrom hash hash values
#' @importFrom rio import
#' @importFrom hash hash
#'
#'   "~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata//NatureBiot2016SuppTable19DoenchRoot.xlsx"
#' source("~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/R/getOfftargetScore2.R")
#' @example
#' mismatch.activity.file <- "~/DropboxUmass/Bioconductor/Trunk/CRISPRseek/inst/extdata/NatureBiot2016SuppTable19DoenchRoot.xlsx"
#'
#' getOfftargetScoreBulge(featureVectors,
#' mismatch.activity.file = mismatch.activity.file)
getOfftargetScoreBulge <-
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
            "NatureBiot2016SuppTable19DoenchRoot.xlsx",
            package = "CRISPRseek")
)
{
    featureVectors <- getOfftargetScore2(featureVectors)
    deletion.activity <- import(mismatch.activity.file, which = 2)
    insertion.activity <- import(mismatch.activity.file, which = 3)
    required.col.ins <- c("Insertion.Type", "Position", "Percent.Active")
    required.col.del <- c("Deletion.Type", "Position", "Percent.Active")
    if (length(intersect(colnames(insertion.activity), required.col.ins)) !=
        length(required.col.ins))
    {
         stop("Please rename the insertion sheet of the activity file column
         to contain at least these 3 column names: Insertion.Type,
              Position, Percent.Active\n")
    }
    if (length(intersect(colnames(deletion.activity), required.col.del)) !=
        length(required.col.del))
    {
        stop("Please rename the deletion sheet of the activity file column
         to contain at least these 3 column names: Deletion.Type,
              Position, Percent.Active\n")
    }
    colnames(insertion.activity)[colnames(insertion.activity) ==
                                     "Insertion.Type"] <- "Type"

    colnames(deletion.activity)[colnames(deletion.activity) ==
                                     "Deletion.Type"] <- "Type"

    #mismatch.activity[mismatch.activity$Mismatch.Type == "rA:dG" &
    #    mismatch.activity$Position == 10,]$Percent.Active
    ##### by default weights is a column vector
    ##### the mismatch activity  is given as pos 20, 19, 18,....1 distance from PAM,    ##### Position named as 1, 2, 3, .....20 though
    ##### and the featureVectors is in the same order now
    ##### so no need to reverse any more. weights = rev(weights)
    fv.geThan1 <- subset(featureVectors, as.numeric(as.character(
        featureVectors$n.insertion)) >= 1)
    fv.lessThan1 <- subset(featureVectors, as.numeric(as.character(
        featureVectors$n.insertion)) < 1)

    featureVectors <- .updateScore(fv.geThan1, fv.lessThan1,
                          col.prefix = "IsInsertion.pos",
                          weights = insertion.activity$Percent.Active,
                          type = insertion.activity$Type,
                          position = insertion.activity$Position)

    # do the same for deletion
    fv.geThan1 <- subset(featureVectors, as.numeric(as.character(
        featureVectors$n.deletion)) >= 1)
    fv.lessThan1 <- subset(featureVectors, as.numeric(as.character(
        featureVectors$n.deletion)) < 1)

    score <- .updateScore(fv.geThan1, fv.lessThan1,
                    col.prefix = "IsDeletion.pos",
                    weights = deletion.activity$Percent.Active,
                    type = deletion.activity$Type,
                    position = deletion.activity$Position)

    score$alignment <- as.character(score$alignment)
    score$score <- round(score$score, 6)
    score  <- score[order(c(score$name,score$score),decreasing = TRUE), ]
    unique(score[!is.na(score$score), ])
}
