#' Default lengthy arguments 
#'
#' This contains a list of long constant values used as defaults in many function.
#' 
#' @format A character string.
#' @examples
#' REpatternFile_default  # Display the default value for REpatternFile.
#' 

#' chromToExclude_default
#' @export
chromToExclude_default = c("chr17_ctg5_hap1", "chr4_ctg9_hap1", "chr6_apd_hap1", "chr6_cox_hap2",  "chr6_dbb_hap3", "chr6_mann_hap4", "chr6_mcf_hap5", "chr6_qbl_hap6", "chr6_ssto_hap7")

#' weights_default
#' @export
weights_default <- c(0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583)

#' gRNA.backbone_default
#' @export
gRNA.backbone_default <- "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU"

#' subPAM.activity_default
#' @export
subPAM.activity_default <- hash(AA = 0, AC = 0, AG = 0.259259259, AT = 0, CA = 0, CC = 0, CG = 0.107142857, CT = 0, GA = 0.069444444, GC = 0.022222222, GG = 1, GT = 0.016129032, TA = 0, TC = 0, TG = 0.038961039, TT = 0)


#' mismatch.activity.file_default
#' @description
#' Default value for mismatch.activity.file (csv), use mismatch.activity.file_default() to show its value
#' @export
mismatch.activity.file_default <- function() {
  return(system.file("extdata", "NatureBiot2016SuppTable19DoenchRoot.csv", package = "CRISPRseek"))
  
}

#' mismatch.activity.file_default_xlsx
#' @description
#' Default value for mismatch.activity.file (xlsx), use mismatch.activity.file_default_xlsx() to show its value.
#' @export
mismatch.activity.file_default_xlsx <- function() {
  return(system.file("extdata", "NatureBiot2016SuppTable19DoenchRoot.xlsx", package = "CRISPRseek"))
}

#' REpatternFile_default
#' @description
#' Default value for REpatternFile, use REpatternFile_default() to show its value.
#' @export
REpatternFile_default <- function() {
  return(system.file("extdata", "NEBenzymes.fa", package = "CRISPRseek"))
}

#' featureWeightMatrixFile_default
#' @description
#' Default value for featureWeightMatrixFile, use featureWeightMatrixFile() to show its value.
#' @export
featureWeightMatrixFile_default <- function() {
  return(system.file("extdata", "DoenchNBT2014.csv", package = "CRISPRseek"))
}

