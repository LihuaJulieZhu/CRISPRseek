#' Output whether the input patterns occurs only once in the sequence
#' 
#' Input a sequence and a list of patterns and determine if the patterns occurs
#' only once in the sequence. Used for determining whether an RE site in gRNA
#' also occurs in the flanking region.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param seq flanking sequence of a gRNA
#' @param patterns patterns as DNAStringSet, such as a list of RE sites
#' @return returns a character vectors containing the uniqueness of each
#' pattern/RE site
#' @note %% ~~further notes~~
#' @author Lihua Julie Zhu
#' @seealso
#' @references %% ~put references to the literature/web site here ~
#' @keywords misc
#' @examples
#' 
#'     seq <- "TGGATTGTATAATCAGCATGGATTTGGAAC"
#'     patterns <- DNAStringSet(c("TGG", "TGGA", "TGGATA", "TTGGAAC", ""))
#'     isPatternUnique(seq, patterns)
#'     isPatternUnique(seq)
#'     isPatternUnique(patterns)
#' @importFrom Biostrings reverseComplement
#' @export
isPatternUnique <- function(seq, patterns)
{
   if (missing(seq) || missing(patterns) || length(seq) == 0 || length(patterns) ==0)
		IsUnique = "NA"
   else
   {
	   IsUnique <- character(length(patterns))
	   seq <- as.character(seq)
	   for (j in 1:length(patterns))
	   {
		   pattern <- patterns[[j]]
		   revpattern <- reverseComplement(pattern)
		   pattern = as.character(pattern)[[1]]
		   if (nchar(pattern) == 0)
		   {
				IsUnique[j] = "NA"
				next
		   }
		   pattern = translatePattern(pattern)
		   revpattern = as.character(revpattern)[[1]]
		   revpattern = translatePattern(revpattern)
		   res.plus <- as.numeric(gregexpr(pattern, seq,
			  perl = TRUE,fixed = FALSE,ignore.case = TRUE)[[1]])
		   n.plus <- length(res.plus[res.plus >0])
		   if (revpattern != pattern)
		   {
		      res.minus <- as.numeric(gregexpr(revpattern, seq,
			     perl = TRUE,fixed = FALSE,ignore.case = TRUE)[[1]])
			  n.minus <- length(res.minus[res.minus >0])
		   }
		   else
			  n.minus <- 0
		   if (n.plus + n.minus > 1)
			   IsUnique[j] <- "No"
		   else if (n.plus + n.minus == 1)
			   IsUnique[j] <- "Yes"
		   else
		       IsUnique[j] <- "NotFound"
	   }
   }
   IsUnique
}
