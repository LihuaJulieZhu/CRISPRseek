#' Predict insertions and deletions induced by CRISPR/Cas9 editing
#'
#' Predict insertions and deletions, and associated reletive frequecies induced by CRISPR/Cas9 editing
#'
#' @details Predict relative indel frequency around target sites of CRISPR/Cas9 system.
#' Currently only Lindel method using logistic regression is implemented in CRISPRseek.
#' Need to download and install the python package Lindel first.
#'
#' git clone https://github.com/shendurelab/Lindel.git 
#'
#' cd Lindel
#'
#' sudo python setup.py install
#'
#' Lindel is compatible with both Python2.7 and Python3.5 or higher.
#'
#' By default, reticulate uses the version of Python found on your PATH (i.e. Sys.which("python")).
#'
#' Set the python path to the output of sys.path (import sys first) in python where Lindel has been setup.
#' For example, use_python('/opt/anaconda3/lib/python3.7')
#'
#' You can also suggest a Conda env to use.
#' For example, use_condaenv("r-cl", conda = "/opt/anaconda3/condabin/conda")
#'
#' Sys.setenv(RETICULATE_PYTHON = "/opt/anaconda3/envs/r-cl/bin/python")
#'
#' \code{}. This function implements the Lindel method
#'
#' @param extendedSequence A vector of DNA sequences of length 60bp. It consists 30bp before the cut site
#' and 30bp after the cut site.

#' @param method the prediction method. default to Lindel. Currently only Lindel method are implemented.
#'
#' @return A list with the same length as the input extendedSequence. 
#'
#' Each list item either contains a data frame with three columns or a warning message.
#'
#' The three columns are the alignment of predicted indel sequence to the original unedited sequence,
#' predicted indel frequency, and the location of the predicted indels.
#' The warning message for Lindel method is as follows.
#'
#' Warning: No PAM sequence is identified. Please check your sequence and try again.
#'
#' @importFrom reticulate py_config py_available py_module_available source_python
#' @export
#'
#' @examples
#' extendedSequence <- c("AAA", "TAACGTTATCAACGCCTATATTAAAGCGACCGTCGGTTGAACTGCGTGGATCAATGCGTC")
#' indelFreq <- predictRelativeFreqIndels(extendedSequence, method = "Lindel") 
#'
#' 
#' @references Wei Chen, Aaron McKenna, Jacob Schreiber et al.,
#' Massively parallel profiling and predictive modeling of the 
#' outcomes of CRISPR/Cas9-mediated double-strand break repair, 
#' Nucleic Acids Research, Volume 47, Issue 15, 05 September 2019, 
#' Pages 7989–8003, https://doi.org/10.1093/nar/gkz487

#' @author Hui Mao and Lihua Julie Zhu

predictRelativeFreqIndels <- function(extendedSequence, method = "Lindel")
{
   if (any(unlist(lapply(extendedSequence, nchar)) != 60))
        warning("Some of the sequences for predicting indel frequencies using Lindel are shorter than 60bp long!")
   py_config()
   if (method == "Lindel")
        e = "Warning: No PAM sequence is identified. Please check your sequence and try again"
   if (! py_available())
	stop("The indel frequency prediction module requires python 2.7, 3.5 or higher to be installed\n")
   pyv <- system2("python", args = "--version", stderr="pythonVersion.txt")
   pyv <- read.table("pythonVersion.txt", sep="", header=FALSE)[1,2]
   pyv <- unlist(strsplit(as.character(pyv), split =".", fixed = TRUE)) 
   if((pyv[1] == 2 && pyv[2] >= 7) || (pyv[1] == 3 && pyv[2] >=5 ))
   {
        unlink("pythonVersion.txt")
        if (! py_module_available("numpy"))
            py_install("numpy")
        if (! py_module_available("scipy"))
            py_install("scipy")
        if( ! py_module_available("pickle")) 
            py_install("pickle")
	if (! py_module_available("re"))
	     py_install("re")
	 if (! py_module_available("json"))
             py_install("json")
        if (! py_module_available("setuptools"))
             py_install("setuptools")
        pyDir <- system.file("extdata/", package = "CRISPRseek")
        oldDir <- getwd()
        if (dir.exists(pyDir))
        {
            setwd(pyDir)
            source_python("RLindelprediction.py")
            tryCatch((
                indelFreq <- lapply(extendedSequence,   function(thisSeq) {
                        res <- predIndelFreq(thisSeq)
                        if (length(res) == 1 && res == e) {
                            return(e)
                        }
                        else {
                            result <- t(data.frame(res))
                            row.names(result) <- NULL
                            colnames(result) <- c("Indels", "Frequency", "Location")
                            return(result)
                         } 
              })), error = function(e) {print(e); })
        }
        else
        {
             stop("The indel frequency prediction scripts are not found in the extdata directory of CRISPRseek, please update CRISPRseek to the newest version first!") 
        }
   }
   else
   {
       stop("Python 2.7, 3.5 or higher is required for predicting relative indel frequency!")
   }
   setwd(oldDir)
   indelFreq 
}
